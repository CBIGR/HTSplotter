import os
import numpy as np


class Savetxt:

    def __init__(self, commonfunctions):
        # self.txt_path, self.date_info, self.fields, self.data, self.elapsed
        # filename, date_info, header, data, enlaspse
        self.header = commonfunctions.fields
        self.date_info = commonfunctions.date_info

        self.data_down = []
        self.data_up =[]

        self.data = np.asarray(commonfunctions.data)
        self.enlapse = np.asarray(commonfunctions.elapsed)
        self.file_path = commonfunctions.txt_path

        # ##atribut from this class

        self.enlapse_complete = None
        self.header_txt_data = None
        self.complete = None
        self.enlap_complete = None
        self.datatype = None

        if len(self.data[0]) == 2:
            self.get_data()
            self.datatype = 1

        self.get_header_format()
        self.get_enlapse()
        self.get_merge_enlapse_header_data()

        self.get_txt_path()

    def get_data(self):
        for i in self.data:
            self.data_down.append(i[0])
            self.data_up.append(i[1])
        self.data_down = np.array(self.data_down)
        self.data_up = np.array(self.data_up)

    def get_txt_path(self):

        # f = open(self.file_path, 'w')
        np.savetxt(self.file_path, self.complete, delimiter='\t', fmt="%s")

        # f.close()

    def get_header_format(self):
        head1 = []
        head2 = []
        head3 = []
        head4 = []
        head5 = []
        head6 = []
        if self.datatype == 1:
            # head1 = [[] for i in range(0, 2)]
            # head2 = [[] for i in range(0, 2)]
            # head3 = [[] for i in range(0, 2)]
            # head4 = [[] for i in range(0, 2)]
            # head5 = [[] for i in range(0, 2)]
            # head6 = [[] for i in range(0, 2)]
            for j in self.header:
                head1.append(j[0])
                head1.append(j[0])
                head2.append(j[1])
                head2.append(j[1])
                head3.append(j[2])
                head3.append(j[2])
                head4.append(j[3])
                head4.append(j[3])
                head5.append(j[4] + "95%CI_down")
                head5.append(j[4] + "95%CI_UP")
                head6.append(j[5])
                head6.append(j[5])
            head1 = np.asarray(head1)
            head2 = np.asarray(head2)
            head3 = np.asarray(head3)
            head4 = np.asarray(head4)
            head5 = np.asarray(head5)
            head6 = np.asarray(head6)

        else:
            for i in self.header:
                head1.append(i[0])
                head2.append(i[1])
                head3.append(i[2])
                head4.append(i[3])
                head5.append(i[4])
                head6.append(i[5])

            head1 = np.asarray(head1)
            head2 = np.asarray(head2)
            head3 = np.asarray(head3)
            head4 = np.asarray(head4)
            head5 = np.asarray(head5)
            head6 = np.asarray(head6)

        head1_2 = np.vstack((head1, head2))
        head1_2_3 = np.vstack((head1_2, head3))
        head1_2_3_4 = np.vstack((head1_2_3, head4))
        head1_2_3_4_5 = np.vstack((head1_2_3_4, head5))
        head1_2_3_4_5_6 = np.vstack((head1_2_3_4_5, head6))

        if self.datatype == 1:
            data_transpodown = self.data_down.T
            data_transpup = self.data_up.T
            data_trans = []
            for k in range(len(self.data_down)):
                data_trans.append(self.data_down[k])
                data_trans.append(self.data_up[k])
            data_transpo = np.asarray(data_trans).T

        else:
            data_transpo = self.data.T
        self.header_txt_data = np.vstack((head1_2_3_4_5_6, data_transpo))

    def get_enlapse(self):

        self.date_info = "_".join(self.date_info)

        head1 = ["Enlapse", "seeding", "condition", "compound", "concentration", self.date_info]

        head1 = np.asarray(head1)
        head1_trans = head1.T

        self.enlap_complete = np.hstack((head1_trans, self.enlapse))

    def get_merge_enlapse_header_data(self):

        self.complete = np.hstack((self.enlap_complete[:, None], self.header_txt_data))


