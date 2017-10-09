def posterior(likelihood, prior):
    def post(x,pks={}):
        return likelihood(x)+prior(x)
    return post