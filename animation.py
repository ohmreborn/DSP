import pandas as pd
import plotly.express as px
df = pd.read_csv("results.csv")
min_range_x = (min(df["x"])-1)*3
max_range_x = (max(df["x"])+1)*3
min_range_y = (min(df["y"])-1)*3
max_range_y = (max(df["y"])+1)*3
range_x = [min_range_x,max_range_x]
range_y = [min_range_y,max_range_y]
fig = px.scatter(df, x="x", y="y", animation_frame="frame", animation_group="star",range_x=range_x,range_y=range_y,color="star")
fig.show()