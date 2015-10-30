#!/usr/bin/env python

"""bm_preproc.py: Boyer-Moore preprocessing."""

__author__ = "Ben Langmead"

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


class TestBoyerMoorePreproc(unittest.TestCase):

    def test_z_1(self):
        s = 'abb'
        #    -00
        z = z_array(s)
        self.assertEqual([3, 0, 0], z)

    def test_z_2(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_z_3(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_n_1(self):
        s = 'abb'
        #    01-
        n = n_array(s)
        self.assertEqual([0, 1, 3], n)

    def test_n_2(self):
        s = 'abracadabra'
        #    1004010100-
        n = n_array(s)
        self.assertEqual([1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11], n)

    def test_n_3(self):
        s = 'abababab'
        #    0204060-
        n = n_array(s)
        self.assertEqual([0, 2, 0, 4, 0, 6, 0, 8], n)

    def test_big_l_prime_1(self):
        s = 'abb'
        #    001
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 2], big_l_prime)

    def test_big_l_prime_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)

    def test_small_l_prime_1(self):
        s = 'abracadabra'
        # N  1004010100-
        # l'           1
        # l'        4
        # l' 44444444111
        small_l_prime = small_l_prime_array(n_array(s))
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

    def test_good_suffix_match_mismatch_1(self):
        p = 'GGTAGGT'
        big_l_prime, big_l, small_l_prime = good_suffix_table(p)
        self.assertEqual([0, 0, 0, 0, 3, 0, 0], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 3, 3, 3], big_l)
        self.assertEqual([7, 3, 3, 3, 3, 0, 0], small_l_prime)
        self.assertEqual(0, good_suffix_mismatch(6, big_l_prime, small_l_prime))
        self.assertEqual(0, good_suffix_mismatch(6, big_l, small_l_prime))
        #  t:      xT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(5, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(5, big_l, small_l_prime))
        #  t:     xGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(4, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(4, big_l, small_l_prime))
        #  t:    xGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(3, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(3, big_l, small_l_prime))
        #  t:   xAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(2, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(2, big_l, small_l_prime))
        #  t:  xTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(1, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(1, big_l, small_l_prime))
        #  t: xGTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(0, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(0, big_l, small_l_prime))

    def test_good_suffix_table_1(self):
        s = 'abb'
        #    001
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 2], big_l_prime)
        self.assertEqual([0, 0, 2], big_l)
        self.assertEqual([3, 0, 0], small_l_prime)

    def test_good_suffix_table_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        # l' -4444444111
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8], big_l)
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

if __name__ == '__main__':
    unittest.main()
