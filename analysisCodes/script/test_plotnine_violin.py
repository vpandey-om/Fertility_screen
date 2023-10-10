import pandas as pd
import numpy as np
import pandas.api.types as pdtypes

from plotnine import *
from plotnine.data import *



np.random.seed(123)
n = 20
mu = (1, 2.3)
sigma = (1, 1.6)

before = np.random.normal(loc=mu[0], scale=sigma[0], size=n)
after = np.random.normal(loc=mu[1], scale=sigma[1], size=n)

df = pd.DataFrame({
    'value': np.hstack([before, after]),
    'when': np.repeat(['before', 'after'], n),
    'id': np.hstack([range(n), range(n)])
})

df['when'] = df['when'].astype(pdtypes.CategoricalDtype(categories=['before', 'after']))

import pdb; pdb.set_trace()
plot=(ggplot(df, aes('when', 'value'))
 + geom_violin(df)
)

plot.save('test.pdf')
