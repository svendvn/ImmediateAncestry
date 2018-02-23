from scipy.spatial.distance import hamming


MOCK_SEQUENCE='eees'+'ssss'+'ssst'+'sstt'+'ettt'+'tttt'+'svsv'+'vvvv'

def mock_likelihood(config, pks={}):
    n=len(config)
    target=MOCK_SEQUENCE[:n]
    return -hamming(list(config),list(target))*1000


if __name__=='__main__':
    print(mock_likelihood('eeee'))