
from textwrap import dedent

from numpy import array


__all__ = [
    'TEST_DNA_STO',
    'TEST_AMINO_STO',
    'TEST_AMINO_NAMES',
    'TEST_STANFEL_NAMES',
    'TEST_Y',
    # TODO: TEST_DNA_X?
    'TEST_AMINO_X',
    'TEST_STANFEL_X'
]


# strip the TEST variables because of the beginning and trailing newlines
TEST_DNA_STO = dedent('''\
    # STOCKHOLM 1.0
    1||A|1       AUGAUUCCCGACUUUAAANNN
    2||A|21      AUGAUUCCCGACUUUAAANNNCAC
    3||A|50      AUGAUUCCCAAANNNCAC
    4||B|0.5     AUGCCCGACUUUAAACAC
    HXB2_env     AUGCCCGACUUUAAACAC
    //''')

TEST_AMINO_STO = dedent('''\
    # STOCKHOLM 1.0
    1||A|1        MIPDFKX-
    2||A|21       MIPDFKXH
    3||A|50       MIP--KXH
    4||B|0.5      .MPDFKH-
    HXB2_env      -MPDFKH-
    //''')

TEST_AMINO_NAMES = ['0aM', '0a[]', 'M1I', 'M1M', 'P2P', 'D3D', 'D3[]', 'F4F', 'F4[]', 'K5K', 'H6H', 'H6X', '6aH', '6a[]']
TEST_STANFEL_NAMES = ['0a[ACGILMPSTV]', '0a[]', 'M1[ACGILMPSTV]', 'P2[ACGILMPSTV]', 'D3[DENQ]', 'D3[]', \
                       'F4[FWY]', 'F4[]', 'K5[HKR]', 'H6[HKR]', 'H6[X]', '6a[HKR]', '6a[]']

TEST_Y = array([1,0,0,1])

TEST_AMINO_X = array([[1,0,1,0,1,1,0,1,0,1,0,1,0,1],
                      [1,0,1,0,1,1,0,1,0,1,0,1,1,0],
                      [1,0,1,0,1,0,1,0,1,1,0,1,1,0],
                      [0,1,0,1,1,1,0,1,0,1,1,0,0,1]])

TEST_STANFEL_X = array([[1,0,1,1,1,0,1,0,1,0,1,0,1],
                        [1,0,1,1,1,0,1,0,1,0,1,1,0],
                        [1,0,1,1,0,1,0,1,1,0,1,1,0],
                        [0,1,1,1,1,0,1,0,1,1,0,0,1]])
