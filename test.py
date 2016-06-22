import progressbar
Nx = 5
Ny = 6
bar = progressbar.ProgressBar(max_value= (Nx + 1)*(Ny + 1))
for i in range(0, Nx + 1):
    for j in range(0, Ny + 1):
        n = (i)*(Nx) + j
        print(i, j , n)
        #bar.update(value = n)