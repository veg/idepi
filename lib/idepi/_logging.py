

__all__ = ['HMMER_LOGGER', '_setup_log']


HMMER_LOGGER = 'FmALfZ4Q4JZ9yVJakdJReEty'


def _setup_log():
    import logging
    h = logging.StreamHandler()
    f = logging.Formatter('%(levelname)s %(asctime)s %(process)d %(funcName)s: %(message)s')
    h.setFormatter(f)
    logging.getLogger(HMMER_LOGGER).addHandler(h)
