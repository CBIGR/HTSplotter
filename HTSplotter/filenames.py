import os


class Filenames:

    def __init__(self, file_name, folderinput, folderresults, information_extracted):
        self.folderinput = folderinput
        self.folderresults = folderresults
        self.information_extracted = information_extracted
        file_name = file_name.split(".")
        self.file_name = file_name[0]
        # output objects
        self.fileitxtnputpath = None
        self.filepngresultspath = None
        self.filepdfresultspath = None
        self.filehdf5resultspath = None
        self.fileictxtresultspath = None
        self.fileicgrtxtresultspath = None
        self.fileblisscorpath = None
        self.filepredictedblisscorpath = None
        self.fileioriginaldatapath = None
        self.fileinhibiteddatapath = None

        self.filegrowthrateresults = None

        self.fileyloewetxtpath = None
        self.fileyloewedifferencetxtpath = None
        self.fileloeweCItxtpath = None

        self.fileyZiptxtpath = None
        self.fileyZipdifferencetxtpath = None
        self.fileZipCItxtpath = None

        self.fileHsascorepath = None
        self.fileHsamaximumpath = None

        self.get_txt_name()
        self.get_informationextracted()
        self.get_error()
        self.get_inhibiteddatatxt_name()
        self.get_normalizedatatxt_name()
        self.get_originaldatatxt_name()
        self.get_ictxt_name()
        self.get_icgrtxt_name()
        self.get_pdf_name()
        self.get_png_name()
        self.get_hdf5_name()

        self.get_blispredictscpretxt_name()
        self.get_blisscpretxt_name()

        self.get_Hsascoretxt_name()
        self.get_Hsamaximumtxt_name()

        self.get_yloewedifference_name()
        self.get_loeweCItxt_name()
        self.get_fileyloewetxt_name()

        self.get_yZipdifference_name()
        self.get_ZipCItxt_name()
        self.get_fileyZiptxt_name()
        self.get_growthratedatatxt_name()

    def get_informationextracted(self):
        termimnation = "_information.txt"
        name = self.file_name + termimnation
        self.information_extractedfile = os.path.join(self.information_extracted, name)

    def get_error(self):
        termimnation = "_Errorfile.txt"
        name = self.file_name + termimnation
        self.errorfile = os.path.join(self.information_extracted, name)

    def get_txt_name(self):
        termimnation = ".txt"
        name = self.file_name + termimnation
        self.fileitxtnputpath = os.path.join(self.folderinput, name)

    def get_ictxt_name(self):
        termimnation = "_IC.txt"
        name = self.file_name + termimnation
        self.fileictxtresultspath = os.path.join(self.folderresults, name)

    def get_icgrtxt_name(self):
        termimnation = "_IC_GrowthRate.txt"
        name = self.file_name + termimnation
        self.fileicgrtxtresultspath = os.path.join(self.folderresults, name)

    def get_blisscpretxt_name(self):
        termimnation = "_Blisscor.txt"
        name = self.file_name + termimnation
        self.fileblisscorpath = os.path.join(self.folderresults, name)

    def get_blispredictscpretxt_name(self):
        termimnation = "_Predicted.txt"
        name = self.file_name + termimnation
        self.filepredictedblisscorpath = os.path.join(self.folderresults, name)

    def get_Hsascoretxt_name(self):
        termimnation = "_HSAscore.txt"
        name = self.file_name + termimnation
        self.fileHsascorepath = os.path.join(self.folderresults, name)

    def get_Hsamaximumtxt_name(self):
        termimnation = "_'HSAmaximum'.txt"
        name = self.file_name + termimnation
        self.fileHsamaximumpath = os.path.join(self.folderresults, name)

    def get_fileyloewetxt_name(self):
        termimnation = "_yloewe.txt"
        name = self.file_name + termimnation
        self.fileyloewetxtpath = os.path.join(self.folderresults, name)

    def get_yloewedifference_name(self):
        termimnation = "_yloewedifference.txt"
        name = self.file_name + termimnation
        self.fileyloewedifferencetxtpath = os.path.join(self.folderresults, name)

    def get_loeweCItxt_name(self):
        termimnation = "_loeweCI.txt"
        name = self.file_name + termimnation
        self.fileloeweCItxtpath = os.path.join(self.folderresults, name)

    def get_fileyZiptxt_name(self):
        termimnation = "_yZip.txt"
        name = self.file_name + termimnation
        self.fileyZiptxtpath = os.path.join(self.folderresults, name)

    def get_yZipdifference_name(self):
        termimnation = "_yZipdifference.txt"
        name = self.file_name + termimnation
        self.fileyZipdifferencetxtpath = os.path.join(self.folderresults, name)

    def get_ZipCItxt_name(self):
        termimnation = "_ZipCI.txt"
        name = self.file_name + termimnation
        self.fileZipCItxtpath = os.path.join(self.folderresults, name)

    def get_originaldatatxt_name(self):
        termimnation = "_Originaldata.txt"
        name = self.file_name + termimnation
        self.fileioriginaldatapath = os.path.join(self.folderresults, name)

    def get_inhibiteddatatxt_name(self):
        termimnation = "_Inhibiteddata.txt"
        name = self.file_name + termimnation
        self.fileinhibiteddatapath = os.path.join(self.folderresults, name)

    def get_normalizedatatxt_name(self):
        termimnation = "_Normalizedatapercentage.txt"
        name = self.file_name + termimnation
        self.filenormalizedatapath = os.path.join(self.folderresults, name)

    def get_png_name(self):
        termimnation = ".png"
        name = self.file_name + termimnation
        self.filepngresultspath = os.path.join(self.folderresults, name)

    def get_pdf_name(self):
        termimnation = ".pdf"
        name = self.file_name + termimnation
        self.filepdfresultspath = os.path.join(self.folderresults, name)

    def get_hdf5_name(self):
        termimnation = ".hdf5"
        name = self.file_name + termimnation
        self.filehdf5resultspath = os.path.join(self.folderresults, name)

    def get_growthratedatatxt_name(self):
        termimnation = "filegrowthrateresults.txt"
        name = self.file_name + termimnation
        self.filegrowthrateresults = os.path.join(self.folderresults, name)
