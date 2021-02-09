__author__ = 'Awab Asher'


from node import Node


class User(Node):

    total_users = 0

    def __init__(self, x_coord, y_coord, coordinates):
        super().__init__(x_coord, y_coord, coordinates)
        self.received_power = None
        self.total_interference = None
        self.bandwidth = None
        self.snr = None
        self.sinr = None
        self.throughput = None
        self.sir = None

        User.total_users += 1

    # getters
    @property
    def get_x_coord(self):
        return self.x_coord

    @property
    def get_y_coord(self):
        return self.y_coord

    @property
    def get_coordinates(self):
        return self.coordinates

    @property
    def get_received_power(self):
        return self.received_power

    @property
    def get_total_interference(self):
        return self.total_interference

    @property
    def get_bandwidth(self):
        return self.bandwidth

    @property
    def get_snr(self):
        return self.snr

    @property
    def get_sinr(self):
        return self.sinr

    @property
    def get_throughput(self):
        return self.throughput

    @property
    def get_sir(self):
        return self.sir

    def __repr__(self):
        return "User(x_coord:{}, y_coord:{}, coordinates:{})".format(self.x_coord,
                                                                     self.y_coord,
                                                                     self.coordinates)

    def __str__(self):
        return "x:{}, " \
               "y:{}, " \
               "coordinates:{}," \
               "received_power:{}," \
               "total_interference:{}," \
               "bandwidth:{}," \
               "snr:{}," \
               "sinr:{}," \
               "throughput:{}," \
               "sir:{}," \
               "".format(self.x_coord,
                         self.y_coord,
                         self.coordinates,
                         self.received_power,
                         self.total_interference,
                         self.bandwidth,
                         self.snr,
                         self.sinr,
                         self.throughput,
                         self.sir)