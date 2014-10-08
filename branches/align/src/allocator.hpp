
#ifndef CCPI_ALIGNED_ALLOCATOR
#define CCPI_ALIGNED_ALLOCATOR

#if LONG_MAX == 2147483647L
typedef unsigned long long ul_int;
#else
typedef unsigned long ul_int;
#endif // LONG_MAX

template <typename T> class aligned_allocator {
public:
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T value_type;

  // must be >= sizeof(pointer) = 8 typically
  // and multiple of which is usually the case with powers of 2
#ifdef NVIDIAGPU
  static const int alignment = 32 * 4; // warp size
#elif defined(__MIC__)
  static const int alignment = 64;
#elif defined(__AVX__)
  static const int alignment = 32;
#else
  static const int alignment = 16;
#endif // x86 alignments

  template<typename U> struct rebind
  { typedef aligned_allocator<U> other; };

  aligned_allocator() { }

  aligned_allocator(const aligned_allocator&) { }

  ~aligned_allocator() { }

  pointer allocate(size_type n, const void *hint = 0)
  {
    std::size_t s = n * sizeof(T) + 2 * alignment;
    void *p = ::operator new(s);
    ul_int offset = ((ul_int)p) % alignment;
    ul_int shift = alignment - offset;
    // if not enough room for a pointer shift into the second block - 2 *
    if (shift < sizeof(pointer))
      shift += alignment;
    void *q = (void *)((ul_int)p + shift);
    pointer r = static_cast<T*>(q);
    // store real pointer just above the one we return
    void **x = (void **)q;
    x = x - 1;
    *x = p;
    return static_cast<T*>(r);
  }

  // p is not permitted to be a null pointer.
  void deallocate(pointer p, size_type)
  {
    void **q = (void **)p;
    q = q - 1;
    void *r = *q;
    ::operator delete(r);
  }

  size_type max_size() const
  {
    // simple way to get 0x8000000000000 without needing L or LL
    return size_type(-1) / sizeof(T);
  }

  void construct(pointer p, const T &val)
  {
    ::new(p) T(val);
  }

  void destroy(pointer p)
  {
    p->~T();
  }
};

// Todo - add an assert #if DEBUG to check that alignment is what it should be?
#if defined(__GNUC__)
#  define assume_aligned(a, T) (T *)__builtin_assume_aligned(a, aligned_allocator<T>::alignment)
#elif defined(__INTEL_COMPILER) || defined(__ICC) || defined(__ICL)
// This doesn't really seem to work. Adding __assume(nv%8==0) etc also failed
// #pragma vector aligned is also supposed to be an alternative
#  define assume_aligned(a,T) a; __assume_aligned(a, aligned_allocator<T>::alignment)
#else
#  define assume_aligned(a,T) a
#endif // compiler choices for forcing alignment
// Todo - need _MSC_VER option

#endif // CCPI_ALIGNED_ALLOCATOR
