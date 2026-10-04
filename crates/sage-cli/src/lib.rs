pub mod input;
pub mod output;
pub mod runner;
pub mod telemetry;

/// Purge memory that the allocator holds freed but committed (mimalloc `mi_collect`,
/// forced). In mimalloc v3 this collects the calling thread's heap and purges the
/// process-wide arenas of the memory already returned to them; pages that other threads
/// still hold in their own heaps are not covered. mimalloc v3 otherwise keeps freed memory
/// committed for up to 1 s (`purge_delay`). `sage` registers this as the database build's
/// release hook ([`sage_core::database::set_release_freed_memory`]), which calls it before
/// the peptide sort and before the fragment index is allocated.
pub fn release_freed_memory() {
    // links the mimalloc library (`libmimalloc-sys`) into every target of this crate,
    // not only into the binary that makes it the global allocator; the Rust binding of
    // `mi_collect` is behind libmimalloc-sys's `extended` feature
    use mimalloc as _;
    extern "C" {
        fn mi_collect(force: bool);
    }
    // SAFETY: `mi_collect` takes no pointers and only returns free memory to the OS.
    unsafe { mi_collect(true) }
}
