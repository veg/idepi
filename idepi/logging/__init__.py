
from __future__ import division, print_function


__all__ = ['IDEPI_LOGGER', 'init_log']


IDEPI_LOGGER = 'FmALfZ4Q4JZ9yVJakdJReEty'


def init_log():
    import logging
    h = logging.StreamHandler()
    f = logging.Formatter('%(levelname)s %(asctime)s %(process)d IDEPI %(funcName)s: %(message)s')
    h.setFormatter(f)
    logging.getLogger(IDEPI_LOGGER).addHandler(h)
