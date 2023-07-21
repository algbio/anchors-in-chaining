import pandas as pd
import matplotlib.pyplot as plt
df = pd.DataFrame({'Animal': ['Falcon', 'Falcon',

                              'Parrot', 'Parrot'],

                   'Max Speed': [380., 370., 24., 26.]})
grouped = df.groupby(['Animal']).mean()
print(grouped)
plt.plot([1,2],grouped['Max Speed'])
plt.show()