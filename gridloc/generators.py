def index_up_down(n_x, n_y):
    y = n_y // 2

    y_sign = -1 if n_y % 2 else 1

    for y_i in range(n_y):

        y += y_sign * y_i
        y_sign = -y_sign

        x_sign = -1 if n_x % 2 else 1
        x = n_x // 2

        for x_i in range(n_x):
            x += x_sign * x_i
            x_sign = -x_sign
            yield (x, y)
