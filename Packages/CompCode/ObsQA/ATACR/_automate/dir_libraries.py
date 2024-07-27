def dir_libraries(CompFolder):
        ATaCR_ML_DataFolder = dict()
        ATaCR_ML_DataFolder['ML_ATaCR_Parent'] = CompFolder + '/ATaCR'
        ATaCR_ML_DataFolder['ML_DataParentFolder'] = CompFolder + '/ATaCR/DATA'
        ATaCR_ML_DataFolder['ML_RawDayData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_day'
        ATaCR_ML_DataFolder['ML_PreProcDayData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_day_preproc'
        ATaCR_ML_DataFolder['ML_RawEventData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_event'
        ATaCR_ML_DataFolder['ML_PreProcEventData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_event_preproc'
        ATaCR_ML_DataFolder['ML_StaSpecAvg'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/AVG_STA'
        ATaCR_ML_DataFolder['ML_CorrectedTraces'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/CORRSEIS'
        ATaCR_ML_DataFolder['ML_b1b2_StaSpectra'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/SPECTRA'
        ATaCR_ML_DataFolder['ML_TransferFunctions'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/TRANSFUN'

        ATaCR_Py_DataFolder = dict()
        ATaCR_Py_DataFolder['Py_DataParentFolder'] = CompFolder + '/ATaCR_Python'
        ATaCR_Py_DataFolder['Py_RawDayData'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/Data'
        # ATaCR_Py_DataFolder['Py_PreProcDayData']
        # ATaCR_Py_DataFolder['Py_RawEventData']
        # ATaCR_Py_DataFolder['Py_PreProcEventData']
        ATaCR_Py_DataFolder['Py_StaSpecAvg'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/AVG_STA'
        ATaCR_Py_DataFolder['Py_CorrectedTraces'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/EVENTS'
        ATaCR_Py_DataFolder['Py_b1b2_StaSpectra'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/SPECTRA'
        ATaCR_Py_DataFolder['Py_TransferFunctions'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/TF_STA'
        ATaCR_Py_DataFolder['Py_Logs'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/Logs'
        return ATaCR_ML_DataFolder,ATaCR_Py_DataFolder
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////