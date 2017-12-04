using DataFrames
mat = [1 2 3 4;3 4 5 6]
df = DataFrame(mat)
writetable("test.csv",df)