# Changelog

## Unreleased

### **Breaking**
- `calcPrefactors` now only returns the prefactors as a Tuple and no longer the mask and number of frequencies. The masking information is now represented by zeros in the prefactors. Previously, masked prefactors where set to one. It is no longer required to output exactly three prefactors
- `mixingFactors` now returns a Matrix of size Fx(N+1) where N is the number of prefactors instead of Fx4. The last column contains the mixing order, while the first N columns contain the mixing factors for the N prefactors