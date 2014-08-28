#ifndef UTILS_H
#define UTILS_H


inline void* safe_malloc(size_t size) {
  void* ptr;
  ptr = malloc(size);
  if (!ptr) {
    fprintf(stderr, "cannot allocate memory");
    exit(7);
  }
  return ptr;
}



#endif // UTILS_H
