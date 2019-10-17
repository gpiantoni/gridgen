
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


def index_spiral(X, Y):
    x = y = 0
    dx = 0
    dy = -1
    for i in range(max(X, Y)**2):
        if (-X / 2 < x <= X / 2) and (-Y / 2 < y <= Y / 2):
            yield (X // 2 - 1 + x, Y // 2 - 1 + y)

        if x == y or (x < 0 and x == -y) or (x > 0 and x == 1 - y):
            dx, dy = -dy, dx
        x, y = x + dx, y + dy
