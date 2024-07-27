def get_ATACR_HPS_Comp(d,atacr_TFs_used='ZP-21', win_length=200):
        RawP = d.Raw[0]['trP'].copy()
        RawZ = d.Raw[0]['trZ'].copy()
        PostATACR = d.Corrected[0][atacr_TFs_used].copy()
        PostATACR.stats.location = 'ATaCR (' + PostATACR.stats.location + ')'
        PostHPS, spectrograms = noisecut(RawZ.copy(), ret_spectrograms=True,win_length=win_length)
        PostBoth, spectrograms = noisecut(PostATACR.copy(), ret_spectrograms=True,win_length=win_length)
        return RawP,RawZ,PostATACR,PostHPS,PostBoth
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
