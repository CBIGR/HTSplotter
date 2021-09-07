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
        self.fileblisscorpath = None
        self.filepredictedblisscorpath = None
        self.fileioriginaldatapath = None
        self.fileinhibiteddatapath = None

        self.get_txt_name()
        self.get_informationextracted()
        self.get_error()
        self.get_inhibiteddatatxt_name()
        self.get_normalizedatatxt_name()
        self.get_originaldatatxt_name()
        self.get_blispredictscpretxt_name()
        self.get_ictxt_name()
        self.get_blisscpretxt_name()
        self.get_pdf_name()
        self.get_png_name()
        self.get_hdf5_name()

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

    def get_blisscpretxt_name(self):
        termimnation = "_Blisscor.txt"
        name = self.file_name + termimnation
        self.fileblisscorpath = os.path.join(self.folderresults, name)

    def get_blispredictscpretxt_name(self):
        termimnation = "_Predicted.txt"
        name = self.file_name + termimnation
        self.filepredictedblisscorpath = os.path.join(self.folderresults, name)

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
