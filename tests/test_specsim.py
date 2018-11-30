import numpy as np
import matplotlib.pyplot as plt
from hyvr.utils import specsim
from hyvr.classes.grid import Grid


def plot_results(twod, covmod):

    x0, y0, z0 = 0., 0., 0.
    dx, dy, dz = 0.1, 0.1, 0.1
    lx, ly, lz = 10., 10., 1.
    grid = Grid(x0, y0, z0, dx, dy, dz, lx, ly, lz, 0)

    var = 1.0
    corl = (1., 2., 0.2)

    if covmod == 'gaussian':
        exponent = 2
    elif covmod == 'exp':
        exponent = 1
    else:
        raise ValueError('No such covmod implemented in test case')

    z = specsim(grid, var, corl, two_dim=twod, covmod=covmod)
    covx = covariogram(z, 0)
    covy = covariogram(z, 1)
    xvec = np.asarray(grid.x[0:len(covx)])
    yvec = np.asarray(grid.y[0:len(covy)])
    covx_theo = var * np.exp(-(xvec/corl[0])**exponent)
    covy_theo = var * np.exp(-(yvec/corl[1])**exponent)
    plt.plot(xvec, covx, 'r-', label=covmod + ' x')
    plt.plot(xvec, covx_theo, 'r--', label=covmod + ' x')
    plt.plot(yvec, covy, 'b-', label=covmod + ' y')
    plt.plot(yvec, covy_theo, 'b--', label=covmod + ' x')
    if not twod:
        covz = covariogram(z, 2)
        zvec = np.asarray(grid.z[0:len(covz)])
        covz_theo = var * np.exp(-(zvec/corl[2])**exponent)
        plt.plot(zvec, covz, 'k-', label=covmod + ' z')
        plt.plot(zvec, covz_theo, 'k--', label=covmod + ' x')
    plt.legend()
    plt.title('var = {:2.4f}'.format(np.var(z)))
    plt.show()

    corrx = np.corrcoef(covx, covx_theo)[0,1]
    corry = np.corrcoef(covy, covy_theo)[0,1]
    assert corrx > 0.5
    assert corry > 0.5
    if not twod:
        corrz = np.corrcoef(covz, covz_theo)[0,1]
        assert corrz > 0.5


def covariogram(z, dim):
    n = z.shape[dim]
    cov = np.zeros(n//2)
    ns = np.zeros(n//2)
    z = z.swapaxes(0, dim)

    for k in range(n//2):
        for i in range(n-k):
            cov[k] += np.sum(z[i,...]*z[i+k,...])
            ns[k] += 1
        ns[k] = z.size/n * (n-k)
    return cov/ns

def test_specsim():
    plot_results(True, 'gaussian')
    plot_results(True, 'exp')
    plot_results(False, 'gaussian')


if __name__ == '__main__':
    test_specsim()
