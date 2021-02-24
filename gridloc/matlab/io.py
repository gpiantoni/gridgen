from mat73 import loadmat as load73


def read_matlab(mat_file):

    out = load73(mat_file)
    if len(out) == 1:
        k = list(out)[0]
        val = out[k]

    else:
        return out

    if k == 'subj_info':
        for field in ('sfile', 'tfile', 'gamma_mean', 'neuralAct'):
            val[field] = (mat_file.parents[1] / val[field])

        val['dims'] = val['dims'].astype(int)

    return val
