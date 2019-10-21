### Create smooth surface

#### Bash
```bash
mris_fill -c -r 1 bert/rh.pial rh_filled.mgz
```

#### Matlab
```matlab
make_outer_surface(...
    'rh_filled.mgz', ...
    15, ...
    'rh.outer', ...
    1)
```

#### Bash
```bash
mris_smooth -nw -n 60 rh.outer rh.smooth
```


