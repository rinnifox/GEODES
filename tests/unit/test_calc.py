import geodes
from pytest import approx


test_structure = 'tests/data/1DB1.pdb'
href = [[127, 142], [149, 152], [226, 246]]


def test_calc_len_of_hel():
    ans = geodes.calc_len_of_hel(test_structure, href)
    correct_ans = [approx(x, rel=0.1) for x in [22.6, 5.6, 30.5]]
    for i, elem in enumerate(correct_ans):
        assert_out_str = (
            "Value {0} at index [{1}] doesn't fall into expected range: {2}"
            .format(
                round(ans[i], 2),
                i,
                elem)
        )
        assert ans[i] == elem, assert_out_str


def test_calc_COM_Calpha_angles():
    ans = geodes.calc_COM_Calpha_angles(test_structure, href)
    correct_ans = [approx(x, rel=0.1) for x in [61.0, 14.2, 120.5]]
    for i, elem in enumerate(correct_ans):
        assert_out_str = (
            "Value {0} at index [{1}] doesn't fall into expected range: {2}"
            .format(
                round(ans[i], 2),
                i,
                elem)
        )
        assert ans[i] == elem, assert_out_str


def test_calc_COM_protein():
    ans = geodes.calc_COM_protein(test_structure)
    correct_ans = [approx(x, rel=0.1) for x in [13.360, 20.776, 42.003]]
    for i, elem in enumerate(correct_ans):
        assert_out_str = (
            "Value {0} at index [{1}] doesn't fall into expected range: {2}"
            .format(
                round(ans[i], 2),
                i,
                elem)
        )
        assert ans[i] == elem, assert_out_str


# def test_calc_angles_between_hel():
#     ans = [calc_angles_between_hel(test_structure, href)]
#     correct_ans = [
#         [approx(49.5, rel=0.1),
#          approx(148.37, rel=0.1)
#         ],
#         approx(125.10, rel=0.1)
#     ]
#     pprint(ans)
#     pprint(correct_ans)
#     assert ans == correct_ans


def test_calc_COM_helix():
    ans = geodes.calc_COM_helix(test_structure, href)
    correct_ans = [[[approx(x, rel=0.1) for x in [18.888, 6.445, 50.056]]],
                   [[approx(x, rel=0.1) for x in [10.596, 10.219, 26.464]]],
                   [[approx(x, rel=0.1) for x in [7.358, 15.986, 35.683]]]]
    for i, elem in enumerate(correct_ans):
        for j, elem2 in enumerate(elem[0]):
            assert_out_str = (
                "Value {0} at index [{1}][0][{2}] doesn't fall into expected range: {3}"
                .format(round(ans[i][0][j], 2),
                        i,
                        j,
                        elem2)
            )
            assert ans[i][0][j] == elem2, assert_out_str


def test_calc_prot_hel_dist():
    ans = geodes.calc_prot_hel_dist(test_structure, href)
    correct_ans = [approx(x, rel=0.1) for x in [17.6, 18.4, 9.9]]
    for i, elem in enumerate(correct_ans):
        assert_out_str = (
            "Value {0} at index [{1}] doesn't fall into expected range: {2}"
            .format(
                round(ans[i], 2),
                i,
                elem)
        )
        assert ans[i] == elem, assert_out_str
