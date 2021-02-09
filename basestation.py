__author__ = 'Awab Asher'

from node import Node


class Basestation(Node):

    total_basestations = 0

    def __init__(self, x_coord, y_coord, coordinates, transmit_power):
        super().__init__(x_coord, y_coord, coordinates)
        self.transmit_power = transmit_power

        Basestation.total_basestations += 1

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
    def get_transmit_power(self):
        return self.transmit_power

    def __repr__(self):
        return "Basestation(x_coord:{}, y_coord:{}, coordinates:{}, transmit_power:{})".format(self.x_coord,
                                                                                               self.y_coord,
                                                                                               self.coordinates,
                                                                                               self.transmit_power)

    def __str__(self):
        return "x:{}, y:{}, coordinates:{}, transmit_power:{}".format(self.x_coord,
                                                                      self.y_coord,
                                                                      self.coordinates,
                                                                      self.transmit_power)