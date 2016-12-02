import pandas as pd

class CCA:
    def __init__(self, x, z, u, v, val_x, val_z):
        self.x = x
        self.z = z
        self.u = u
        self.v = v
        self.val_x = val_x
        self.val_z = val_z

    def num_nonzero_components(self):
        # count # of nonzero components in u and v
        pass

    def project_x(self):
        return self.x.dot(self.u)

    def project_val_x(self):
        return self.val_x.dot(self.u)

    def project_z(self):
        return self.z.dot(self.v)

    def project_val_z(self):
        return self.val_z.dot(self.v)

    def plot_projections(self):
        # scatter plot of x.dot(u) vs z.dot(v)
        # one color for train, one color for validation
        pass

    @staticmethod
    def correlation(x, z):
        # return correlation of projections for x.dot(u) and z.dot(v)
        pass

    def summarise(self):
        summary = {}
        summary['train correlation'] = self.correlation(self.x, self.z)
        summary['validation correlation'] = \
            self.correlation(self.val_x, self.val_z)
        summary['# nonzero u weights'] = None
        summary['# nonzero v weights'] = None

        summary = {k:[v] for k, v in summary}

        return pd.DataFrame()
