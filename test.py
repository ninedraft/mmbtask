import progressbar

bar = progressbar.ProgressBar(max_value= 5*6)
for i in range(0, 5):
    for j in range(0, 6):
        bar.update((i + 1)*5 + j)