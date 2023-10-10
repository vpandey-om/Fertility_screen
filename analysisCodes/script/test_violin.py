from plotnine import ggplot, aes, geom_point, facet_wrap, geom_text
import pandas as pd

# Create a sample DataFrame
df = pd.DataFrame({'x': [1, 2, 3, 1, 2, 3],
                   'y': [4, 5, 6, 7, 8, 9],
                   'facet_var': ['A', 'A', 'A', 'B', 'B', 'B'],
                   'label': ['Text A1', 'Text A2', 'Text A3', 'Text B1', 'Text B2', 'Text B3']})

# Create the plot using Plotnine
plot =ggplot(df, aes(x='x', y='y')) + geom_point()

# Add text annotations to each facet
plot += facet_wrap('~ facet_var', scales='free')

# Custom annotations for specific facets
custom_annotations = pd.DataFrame({
    'facet_var': ['A'],
    'x': [1.5],
    'y': [8],
    'label': ['Custom A']
})

# Add custom annotations to specific facets
plot += geom_text( aes(x='x', y='y', label='label'), data=custom_annotations,size=10, color='blue')

# Display the plot
print(plot)



plot.save('test.pdf')
