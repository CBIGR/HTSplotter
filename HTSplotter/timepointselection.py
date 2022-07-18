import math


class Timepointselection:
    def __init__(self, elapse):
        self.elapse = elapse

        self.time_selected = []
        self.time_position = []

        self.time_point_position()

    def time_point_position(self):

        if len(self.elapse) > 1:
            main_time = math.ceil(24 / self.elapse[1])
            # scans_per_day = 24 / x_main_times
            if main_time == 1:  # each 24h interval time points
                x_main_times = len(self.elapse) - 1
            else:
                x_main_times = round(len(self.elapse) / (main_time + 1))
            count = 1
            for i in range(x_main_times):
                try:
                    self.time_selected.append(self.elapse[count * main_time])
                    self.time_position.append(count * main_time)
                except IndexError:
                    self.time_selected.append(self.elapse[-1])
                    self.time_position.append(len(self.elapse) - 1)
                    pass

                count += 1
        else:
            self.time_selected.append(self.elapse[0])
            self.time_position.append(0)
