import scipy.stats
import numpy as np

def lin_fit(x,y):
    result = scipy.stats.linregress(x, y)
    tinv = lambda p, df: abs(scipy.stats.t.ppf(p/2, df))
    ts = tinv(0.05, len(x)-2)
    slpe_ci=np.array([0.0,0.0]);intrcpt_ci=np.array([0.0,0.0])
    slpe_ci[0] = result.slope - result.stderr * ts
    slpe_ci[1] = result.slope + result.stderr * ts
    intrcpt_ci[0] = result.intercept - result.intercept_stderr * ts
    intrcpt_ci[1] = result.intercept + result.intercept_stderr * ts
    R2=result.rvalue**2
    t = np.sqrt(R2 * (len(x) - 2)) / np.sqrt(1 - R2)
    t1=scipy.stats.t.ppf(0.975, len(x) - 2)
    t2 = scipy.stats.t.ppf(0.995, len(x) - 2)
    t3 = scipy.stats.t.ppf(0.9995, len(x) - 2)

    if t < t1:
        signif = 0;signif_label=''
    elif t < t2:
        signif = 1;signif_label='*'
    elif t < t3:
        signif = 2;signif_label='**'
    elif t >= t3:
        signif = 3;signif_label='***'
    return result,slpe_ci,intrcpt_ci,signif,signif_label

if __name__ == "__main__":
    x = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    y = np.array([2, 1, 4, 5, 8, 12, 18, 25, 96, 48])
    (result,slpe_ci,intrcpt_ci,signif, signif_label) = lin_fit(x, y)
    print(result)
    print('significativity of fit is %s (%s)' % signif,signif_label)