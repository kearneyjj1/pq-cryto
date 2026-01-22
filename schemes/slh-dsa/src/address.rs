//! Address (ADRS) structure for SLH-DSA domain separation.
//!
//! The ADRS structure is a 32-byte value used to ensure domain separation
//! for all hash function calls in SLH-DSA. This prevents related-key attacks
//! and ensures that different parts of the scheme use different hash instances.
//!
//! # Structure (FIPS 205 Section 4.2)
//!
//! ```text
//! Bytes 0-3:   Layer address
//! Bytes 4-15:  Tree address (96 bits)
//! Bytes 16-19: Type
//! Bytes 20-31: Type-specific data (depends on address type)
//! ```

/// Address types for SLH-DSA (FIPS 205 Section 4.2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u32)]
pub enum AddressType {
    /// WOTS+ hash address (for chain computations)
    WotsHash = 0,
    /// WOTS+ public key compression address
    WotsPk = 1,
    /// Tree hash address (for internal Merkle nodes)
    Tree = 2,
    /// FORS tree leaf generation address
    ForsTree = 3,
    /// FORS roots compression address
    ForsPk = 4,
    /// WOTS+ key generation address
    WotsPrf = 5,
    /// FORS key generation address
    ForsPrf = 6,
}

/// ADRS structure (32 bytes) for domain separation in SLH-DSA.
///
/// This structure ensures that every hash call in SLH-DSA uses a unique
/// domain separator, preventing any unintended relationships between
/// different hash evaluations.
#[derive(Clone, Copy, Debug)]
pub struct Address {
    /// Raw 32-byte address data
    data: [u8; 32],
}

impl Default for Address {
    fn default() -> Self {
        Self::new()
    }
}

impl Address {
    /// Creates a new zeroed address.
    pub fn new() -> Self {
        Address { data: [0u8; 32] }
    }

    /// Sets the layer address (bytes 0-3).
    ///
    /// The layer indicates which level in the hypertree we're operating on.
    /// Layer 0 is the bottom, layer d-1 is the top.
    #[inline]
    pub fn set_layer(&mut self, layer: u32) {
        self.data[0..4].copy_from_slice(&layer.to_be_bytes());
    }

    /// Gets the layer address.
    #[inline]
    pub fn layer(&self) -> u32 {
        u32::from_be_bytes([self.data[0], self.data[1], self.data[2], self.data[3]])
    }

    /// Sets the tree address (bytes 4-15, 96 bits).
    ///
    /// This identifies which tree within the layer we're operating on.
    /// Only the lower 96 bits of the u128 are used.
    #[inline]
    pub fn set_tree(&mut self, tree: u64) {
        // Tree address is 12 bytes (96 bits), but we use u64 for simplicity
        // Store in bytes 4-11, leaving bytes 12-15 as zeros
        self.data[4..12].copy_from_slice(&tree.to_be_bytes());
        self.data[12..16].fill(0);
    }

    /// Gets the tree address as u64.
    #[inline]
    pub fn tree(&self) -> u64 {
        u64::from_be_bytes([
            self.data[4],
            self.data[5],
            self.data[6],
            self.data[7],
            self.data[8],
            self.data[9],
            self.data[10],
            self.data[11],
        ])
    }

    /// Sets the address type (bytes 16-19).
    #[inline]
    pub fn set_type(&mut self, addr_type: AddressType) {
        let type_val = addr_type as u32;
        self.data[16..20].copy_from_slice(&type_val.to_be_bytes());
        // Clear type-specific data when type changes
        self.data[20..32].fill(0);
    }

    /// Gets the address type.
    #[inline]
    pub fn addr_type(&self) -> u32 {
        u32::from_be_bytes([self.data[16], self.data[17], self.data[18], self.data[19]])
    }

    // ========================================================================
    // WOTS+ specific fields (for WotsHash, WotsPk, WotsPrf types)
    // ========================================================================

    /// Sets the key pair address (bytes 20-23).
    ///
    /// Identifies which WOTS+ key pair within the XMSS tree.
    #[inline]
    pub fn set_keypair(&mut self, keypair: u32) {
        self.data[20..24].copy_from_slice(&keypair.to_be_bytes());
    }

    /// Gets the key pair address.
    #[inline]
    pub fn keypair(&self) -> u32 {
        u32::from_be_bytes([self.data[20], self.data[21], self.data[22], self.data[23]])
    }

    /// Sets the chain address (bytes 24-27).
    ///
    /// Identifies which chain within the WOTS+ signature.
    #[inline]
    pub fn set_chain(&mut self, chain: u32) {
        self.data[24..28].copy_from_slice(&chain.to_be_bytes());
    }

    /// Gets the chain address.
    #[inline]
    pub fn chain(&self) -> u32 {
        u32::from_be_bytes([self.data[24], self.data[25], self.data[26], self.data[27]])
    }

    /// Sets the hash address (bytes 28-31).
    ///
    /// Identifies which hash within the chain.
    #[inline]
    pub fn set_hash(&mut self, hash: u32) {
        self.data[28..32].copy_from_slice(&hash.to_be_bytes());
    }

    /// Gets the hash address.
    #[inline]
    pub fn hash(&self) -> u32 {
        u32::from_be_bytes([self.data[28], self.data[29], self.data[30], self.data[31]])
    }

    // ========================================================================
    // Tree/FORS specific fields (for Tree, ForsTree, ForsPk, ForsPrf types)
    // ========================================================================

    /// Sets the tree height (bytes 24-27).
    ///
    /// Used for Tree and ForsTree types to indicate the height in the tree.
    #[inline]
    pub fn set_tree_height(&mut self, height: u32) {
        self.data[24..28].copy_from_slice(&height.to_be_bytes());
    }

    /// Gets the tree height.
    #[inline]
    pub fn tree_height(&self) -> u32 {
        u32::from_be_bytes([self.data[24], self.data[25], self.data[26], self.data[27]])
    }

    /// Sets the tree index (bytes 28-31).
    ///
    /// Identifies which node at the given height.
    #[inline]
    pub fn set_tree_index(&mut self, index: u32) {
        self.data[28..32].copy_from_slice(&index.to_be_bytes());
    }

    /// Gets the tree index.
    #[inline]
    pub fn tree_index(&self) -> u32 {
        u32::from_be_bytes([self.data[28], self.data[29], self.data[30], self.data[31]])
    }

    /// Returns the address as a 32-byte array for use in hash functions.
    #[inline]
    pub fn to_bytes(&self) -> [u8; 32] {
        self.data
    }

    /// Returns a reference to the internal byte array.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.data
    }

    /// Creates a copy of this address for use in nested operations.
    pub fn copy(&self) -> Self {
        Address { data: self.data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_address_is_zeroed() {
        let addr = Address::new();
        assert_eq!(addr.data, [0u8; 32]);
    }

    #[test]
    fn test_set_layer() {
        let mut addr = Address::new();
        addr.set_layer(7);
        assert_eq!(addr.layer(), 7);
        assert_eq!(&addr.data[0..4], &[0, 0, 0, 7]);
    }

    #[test]
    fn test_set_tree() {
        let mut addr = Address::new();
        addr.set_tree(0x123456789ABCDEF0);
        assert_eq!(addr.tree(), 0x123456789ABCDEF0);
    }

    #[test]
    fn test_set_type_clears_data() {
        let mut addr = Address::new();
        addr.set_keypair(42);
        addr.set_chain(10);
        addr.set_hash(5);

        // Setting type should clear type-specific data
        addr.set_type(AddressType::WotsHash);
        assert_eq!(addr.keypair(), 0);
        assert_eq!(addr.chain(), 0);
        assert_eq!(addr.hash(), 0);
    }

    #[test]
    fn test_wots_fields() {
        let mut addr = Address::new();
        addr.set_type(AddressType::WotsHash);
        addr.set_keypair(100);
        addr.set_chain(15);
        addr.set_hash(7);

        assert_eq!(addr.keypair(), 100);
        assert_eq!(addr.chain(), 15);
        assert_eq!(addr.hash(), 7);
    }

    #[test]
    fn test_tree_fields() {
        let mut addr = Address::new();
        addr.set_type(AddressType::Tree);
        addr.set_keypair(50); // Using keypair field for tree operations
        addr.set_tree_height(3);
        addr.set_tree_index(12);

        assert_eq!(addr.keypair(), 50);
        assert_eq!(addr.tree_height(), 3);
        assert_eq!(addr.tree_index(), 12);
    }

    #[test]
    fn test_address_types() {
        assert_eq!(AddressType::WotsHash as u32, 0);
        assert_eq!(AddressType::WotsPk as u32, 1);
        assert_eq!(AddressType::Tree as u32, 2);
        assert_eq!(AddressType::ForsTree as u32, 3);
        assert_eq!(AddressType::ForsPk as u32, 4);
        assert_eq!(AddressType::WotsPrf as u32, 5);
        assert_eq!(AddressType::ForsPrf as u32, 6);
    }

    #[test]
    fn test_to_bytes() {
        let mut addr = Address::new();
        addr.set_layer(1);
        addr.set_tree(2);
        addr.set_type(AddressType::WotsHash);
        addr.set_keypair(3);
        addr.set_chain(4);
        addr.set_hash(5);

        let bytes = addr.to_bytes();
        assert_eq!(bytes.len(), 32);

        // Verify layer
        assert_eq!(&bytes[0..4], &[0, 0, 0, 1]);
        // Verify type
        assert_eq!(&bytes[16..20], &[0, 0, 0, 0]); // WotsHash = 0
    }

    #[test]
    fn test_copy() {
        let mut addr = Address::new();
        addr.set_layer(5);
        addr.set_keypair(10);

        let addr_copy = addr.copy();
        assert_eq!(addr.layer(), addr_copy.layer());
        assert_eq!(addr.keypair(), addr_copy.keypair());
    }
}
