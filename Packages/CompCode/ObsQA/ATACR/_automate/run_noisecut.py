def run_noisecut(tr,win_length=200,width=None):
        hps, (S_full, S_background, S_hps, frequencies, times) = noisecut(tr, ret_spectrograms=True,win_length=win_length,width=width)
        return hps, (S_full, S_background, S_hps, frequencies, times)