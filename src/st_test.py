"""Test suffix trees code."""

import collections
from typing import (Callable, Iterator)

from test_helpers import (
    _Test, check_equal_matches, check_matches, check_sorted,
    collect_tests, fibonacci_string, pick_random_patterns,
    pick_random_patterns_len, pick_random_prefix,
    pick_random_suffix, random_string
)

from alphabet import (Alphabet)
from st import (
    Inner, Leaf, SuffixTree,
    mccreight_st_construction, naive_st_construction
)

STConstructor = Callable[[str], SuffixTree]


ALGOS: list[STConstructor] = [
    naive_st_construction,
    mccreight_st_construction,
]


def strip_algo_name(name: str) -> str:
    """Get the algorithm name from the function name."""
    return name.split('_')[0]


def test_contains() -> None:
    """Check a suffix tree's contain method."""
    st = naive_st_construction("mississippi")
    assert "iss" in st
    assert "sss" not in st
    assert "ip" in st
    assert "x" not in st

    st = mccreight_st_construction("mississippi")
    assert "iss" in st
    assert "sss" not in st
    assert "ip" in st
    assert "x" not in st


def check_st_sorted(algo: STConstructor) -> _Test:
    """Check that suffixes are sorted."""
    def test(_: object) -> None:
        for _ in range(10):
            x = random_string(20, alpha="abc")
            # using the leaf iterator
            check_sorted(x, list(algo(x).root))
        for fib in range(5, 10):
            x = fibonacci_string(fib)
            check_sorted(x, list(algo(x).root))
        for n in range(5, 50):
            x = 'a' * n
            check_sorted(x, list(algo(x).root))
    return test


TestSorted = collect_tests(
    (strip_algo_name(algo.__name__), check_st_sorted(algo))
    for algo in ALGOS
)


def test_node_comparison() -> None:
    """Test that we can compare nodes."""
    x, _ = Alphabet.mapped_subseq("foo")
    y, _ = Alphabet.mapped_subseq("bar")
    assert Leaf(0, x) == Leaf(0, x)
    assert Leaf(0, x) != Leaf(1, x)
    assert Leaf(0, x) != Leaf(0, y)

    assert Inner(x) == Inner(x)
    assert Inner(x) != Inner(y)

    inner1 = Inner(x)
    inner2 = Inner(x)
    assert inner1 == inner2
    inner1.add_children(Leaf(0, y))
    assert inner1 != inner2
    inner2.add_children(Leaf(0, y))
    assert inner1 == inner2


def check_equal_mccreight(algo: STConstructor) -> _Test:
    """Check that constructions are equal to McCreight's tree."""
    def test(_: object) -> None:
        for _ in range(10):
            x = random_string(20, alpha="abc")
            assert mccreight_st_construction(x) == algo(x)
        for fib in range(5, 10):
            x = fibonacci_string(fib)
            assert mccreight_st_construction(x) == algo(x)
        for n in range(5, 50):
            x = 'a' * n
            assert mccreight_st_construction(x) == algo(x)
    return test


TestEqualMcCreight = collect_tests(
    (strip_algo_name(algo.__name__), check_equal_mccreight(algo))
    for algo in ALGOS
)


def check_occurrences(algo: STConstructor) -> _Test:
    """Check searching with suffix trees."""
    def st_search(x: str, p: str) -> Iterator[int]:
        return algo(x).search(p)

    def test(_: object) -> None:
        for _ in range(10):
            x = random_string(50, alpha="abcd")
            for p in pick_random_patterns(x, 5):
                check_matches(x, p, st_search(x, p))
            for p in pick_random_patterns_len(x, 5, 3):
                check_matches(x, p, st_search(x, p))
            for p in pick_random_prefix(x, 5):
                check_matches(x, p, st_search(x, p))
            for p in pick_random_suffix(x, 5):
                check_matches(x, p, st_search(x, p))
    return test


TestMatches = collect_tests(
    (strip_algo_name(algo.__name__), check_occurrences(algo))
    for algo in ALGOS
)


def bmh(x: str, p: str) -> Iterator[int]:
    """Run the Boyer-Moore-Horspool algorithm."""
    # Can't handle empty strings directly
    if not p:
        yield from range(len(x) + 1)
        return

    jump: dict[str, int] = collections.defaultdict(lambda: len(p))
    for j, a in enumerate(p[:-1]):  # skip last index!
        jump[a] = len(p) - j - 1

    i, j = 0, 0
    while i < len(x) - len(p) + 1:
        for j in reversed(range(len(p))):
            if x[i + j] != p[j]:
                break
        else:
            yield i

        i += jump[x[i + len(p) - 1]]


def check_against_bmh(algo: STConstructor) -> _Test:
    """Check that suffix trees find the same matches as BMH."""
    def st_search(x: str, p: str) -> Iterator[int]:
        print(algo(x).search(p))
        print(list(algo(x).search(p)))
        return algo(x).search(p)

    def test(_: object) -> None:
        for _ in range(10):
            x = random_string(50, alpha="abcd")
            for p in pick_random_patterns(x, 5):
                check_equal_matches(x, p, bmh, st_search)
            for p in pick_random_patterns_len(x, 5, 3):
                check_equal_matches(x, p, bmh, st_search)
            for p in pick_random_prefix(x, 5):
                check_equal_matches(x, p, bmh, st_search)
            for p in pick_random_suffix(x, 5):
                check_equal_matches(x, p, bmh, st_search)
    return test


TestAgainstBMH = collect_tests(
    (strip_algo_name(algo.__name__), check_against_bmh(algo))
    for algo in ALGOS
)
