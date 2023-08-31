import numpy as np
from scipy import linalg, stats
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
import os, sys, math

import ray


from .utility import CPM



@ray.remote(num_cpus=2)
class IRMAD:


    def __init__(self, master, slave, output, filename, out_fmt='GTiff', penalization=0.001, logger=None):

        gdal.AllRegister()

        self.logger = logger
        if output:
            os.chdir(output)

        master_ds = gdal.Open(master, GA_ReadOnly)
        slave_ds = gdal.Open(slave, GA_ReadOnly)
        if not master_ds and slave_ds:
            pass
        else:

            self.master_cols, self.master_rows, self.master_bands, self.master_dataset = self._get_raster_metadata(master)
            self.master_pos = range(1, self.master_bands + 1)

            self.slave_cols, self.slave_rows, self.slave_bands, self.slave_dataset = self._get_raster_metadata(slave)
            self.slave_pos = range(1, self.slave_bands + 1)

            # penalization
            self.lam = penalization

            self.fmt = out_fmt
            self.outfile = filename

            # match dimensions
            self.cols, self.rows, self.bands = self._get_dims()

        # else:
        #     sys.stderr.write("Please import two images as master and slave to perform iMad processing")
        #     sys.exit(1)

    def _get_raster_metadata(self, image):

        inDataset = gdal.Open(image, GA_ReadOnly)
        cols = inDataset.RasterXSize
        rows = inDataset.RasterYSize
        bands = inDataset.RasterCount

        return (cols, rows, bands, inDataset)

    def _get_dims(self):
        if (self.master_cols != self.slave_cols) or (self.master_rows != self.slave_rows)\
                or (self.master_bands != self.slave_bands) or (len(self.master_pos) != len(self.slave_pos)):
            sys.stderr.write("Size mismatch")
            sys.exit(1)
        else:
            return self.master_cols, self.master_rows, self.master_bands

    def _choldc(self, A):
        # Cholesky-Banachiewicz algorithm,
        # A is a numpy matrix
        L = A - A
        for i in range(len(L)):
            for j in range(i):
                sm = 0.0
                for k in range(j):
                    sm += L[i, k] * L[j, k]
                L[i, j] = (A[i, j] - sm) / L[j, j]
            sm = 0.0
            for k in range(i):
                sm += L[i, k] * L[i, k]
            L[i, i] = math.sqrt(A[i, i] - sm)
        return L

    def _geneiv(self, A, B):
        # solves A*x = lambda*B*x for numpy matrices A and B,
        # returns eigenvectors in columns
        Li = np.linalg.inv(self._choldc(B))
        C = Li * A * (Li.transpose())
        C = np.asmatrix((C + C.transpose()) * 0.5, np.float32)
        eivs, V = np.linalg.eig(C)
        return eivs, Li.transpose() * V

    def MAD_iteration(self):

        if not os.path.exists(self.outfile):
            cpm = CPM(2 * self.bands)
            delta = 1.0
            oldrho = np.zeros(self.bands)
            itr = 0
            tile = np.zeros((self.cols, 2 * self.bands))
            sigMADs = 0
            means1 = 0
            means2 = 0
            A = 0
            B = 0
            rasterBands1 = []
            rasterBands2 = []
            for b in self.master_pos:
                rasterBands1.append(self.master_dataset.GetRasterBand(b))
            for b in self.master_pos:
                rasterBands2.append(self.slave_dataset.GetRasterBand(b))

            while (delta > 0.001) and (itr < 100):
                # spectral tiling for statistics
                for row in range(self.rows):
                    for k in range(self.bands):
                        tile[:, k] = rasterBands1[k].ReadAsArray(0, 0 + row, self.cols, 1)
                        tile[:, self.bands + k] = rasterBands2[k].ReadAsArray(0, 0 + row, self.cols, 1)

                    # eliminate no-data pixels (assuming all zeroes)
                    tst1 = np.sum(tile[:, 0:self.bands], axis=1)
                    tst2 = np.sum(tile[:, self.bands::], axis=1)
                    idx1 = set(np.where((tst1 > 0))[0])
                    idx2 = set(np.where((tst2 > 0))[0])
                    idx = list(idx1.intersection(idx2))
                    if itr > 0:
                        mads = np.asarray((tile[:, 0:self.bands] - means1) * A - (tile[:, self.bands::] - means2) * B)
                        chisqr = np.sum((mads / sigMADs) ** 2, axis=1)
                        wts = 1 - stats.chi2.cdf(chisqr, [self.bands])
                        cpm.update(tile[idx, :], wts[idx])
                    else:
                        cpm.update(tile[idx, :])

                # weighted covariance matrices and means
                S = cpm.covariance()
                means = cpm.means()

                # reset prov means object
                CPM(2 * self.bands)
                s11 = S[0:self.bands, 0:self.bands]
                s11 = (1 - self.lam) * s11 + self.lam * np.eye(self.bands)
                s22 = S[self.bands:, self.bands:]
                s22 = (1 - self.lam) * s22 + self.lam * np.eye(self.bands)
                s12 = S[0:self.bands, self.bands:]
                s21 = S[self.bands:, 0:self.bands]
                c1 = s12 * linalg.inv(s22) * s21
                b1 = s11
                c2 = s21 * linalg.inv(s11) * s12
                b2 = s22

                # solution of generalized eigen problems
                if self.bands > 1:
                    mu2a, A = self._geneiv(c1, b1)
                    mu2b, B = self._geneiv(c2, b2)

                    # sort a
                    idx = np.argsort(mu2a)
                    A = A[:, idx]

                    # sort b
                    idx = np.argsort(mu2b)
                    B = B[:, idx]
                    mu2 = mu2b[idx]

                else:
                    mu2 = c1 / b1
                    A = 1 / np.sqrt(b1)
                    B = 1 / np.sqrt(b2)

                # canonical correlations
                mu = np.sqrt(mu2)
                a2 = np.diag(A.T * A)
                b2 = np.diag(B.T * B)
                sigma = np.sqrt((2 - self.lam * (a2 + b2)) / (1 - self.lam) - 2 * mu)
                rho = mu * (1 - self.lam) / np.sqrt((1 - self.lam * a2) * (1 - self.lam * b2))

                # stopping criterion
                delta = max(abs(rho - oldrho))
                print(delta, rho)
                oldrho = rho

                # tile the sigmas and means
                sigMADs = np.tile(sigma, (self.cols, 1))
                means1 = np.tile(means[0:self.bands], (self.cols, 1))
                means2 = np.tile(means[self.bands::], (self.cols, 1))

                # ensure sum of positive correlations between X and U is positive
                D = np.diag(1 / np.sqrt(np.diag(s11)))
                s = np.ravel(np.sum(D * s11 * A, axis=0))
                A = A * np.diag(s / np.abs(s))

                # ensure positive correlation between each pair of canonical variates
                cov = np.diag(A.T * s12 * B)
                B = B * np.diag(cov / np.abs(cov))
                itr += 1


            # write results to disk
            driver = gdal.GetDriverByName(self.fmt)
            outDataset = driver.Create(self.outfile, self.cols, self.rows, self.bands + 1, GDT_Float32, options=['COMPRESS=LZW'])
            projection = self.master_dataset.GetProjection()
            geotransform = self.master_dataset.GetGeoTransform()

            if geotransform is not None:
                gt = list(geotransform)
                gt[0] = gt[0] + 0 * gt[1]
                gt[3] = gt[3] + 0 * gt[5]
                outDataset.SetGeoTransform(tuple(gt))
            if projection is not None:
                outDataset.SetProjection(projection)

            outBands = []
            for k in range(self.bands + 1):
                outBands.append(outDataset.GetRasterBand(k + 1))

            for row in range(self.rows):
                for k in range(self.bands):
                    tile[:, k] = rasterBands1[k].ReadAsArray(0, 0 + row, self.cols, 1)
                    tile[:, self.bands + k] = rasterBands2[k].ReadAsArray(0, 0 + row, self.cols, 1)
                mads = np.asarray((tile[:, 0:self.bands] - means1) * A - (tile[:, self.bands::] - means2) * B)
                chisqr = np.sum((mads / sigMADs) ** 2, axis=1)

                for k in range(self.bands):
                    outBands[k].WriteArray(np.reshape(mads[:, k], (1, self.cols)), 0, row)
                outBands[self.bands].WriteArray(np.reshape(chisqr, (1, self.cols)), 0, row)

            for outBand in outBands:
                outBand.FlushCache()
            outDataset = None
            inDataset1 = None
            inDataset2 = None
            self.logger.info('result written to: ' + self.outfile)
            self.logger.info('-----------------done---------------------')
