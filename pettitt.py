import numpy as np

def pettitt_test(X):
    """
    Pettitt test calculated following Pettitt (1979): https://www.jstor.org/stable/2346729?seq=4#metadata_info_tab_contents
    """

    T = len(X)
    U = []
    for t in range(T): # t is used to split X into two subseries
        X_stack = np.zeros((t, len(X[t:]) + 1), dtype=int)
        X_stack[:,0] = X[:t] # first column is each element of the first subseries
        X_stack[:,1:] = X[t:] # all rows after the first element are the second subseries
        U.append(np.sign(X_stack[:,0] - X_stack[:,1:].transpose()).sum()) # sign test between each element of the first subseries and all elements of the second subseries, summed.

    tau = np.argmax(np.abs(U)) # location of change (first data point of second sub-series)
    K = np.max(np.abs(U))
    p = 2 * np.exp(-6 * K**2 / (T**3 + T**2))
        
    return (tau, p)