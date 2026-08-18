# cghash
Contains implementations of cryptographic hash functions into class groups of binary quadratic forms. 

## hash algorithms
The basic implementation **pr** is based on [Wesolowski's proposal](https://doi.org/10.1007/s00145-020-09364-x). It does rejection sampling of the first coefficient until it passes a Miller-Rabin test with n=25. The variation **pp** is derived from [Schoenmakers' and Segers' encoding](https://doi.org/10.3390/cryptography7040056), where we immediately compute the second coefficient without testing for primality. Lastly, **pp-mr2** is a hybrid of the two, where a base-2 Miller-Rabin test is done before computing the second coefficient.

Computing the remainder of Δ modulo the candidate coefficient at the right moment is an important optimisation for the **pp** and **pp-mr2** variants. This technique is denoted with **-dma**.

Finally, the prefix **clr-** indicates hashing to multiple forms and composing, as suggested by [Chalkias, Lindstrøm and Røy](https://eprint.iacr.org/2024/295).

