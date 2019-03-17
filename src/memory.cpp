//
// Created by xiamr on 3/17/19.
//

#include <exception>
#include <tbb/scalable_allocator.h>

void *operator new(std::size_t size) noexcept(false) {
    if (size == 0) size = 1;
    if (void *ptr = scalable_malloc(size))
        return ptr;
    throw std::bad_alloc();
}

void *operator new[](std::size_t size) noexcept(false) {
    return operator new(size);
}

void *operator new(std::size_t size, const std::nothrow_t &) noexcept {
    if (size == 0) size = 1;
    if (void *ptr = scalable_malloc(size))
        return ptr;
    return nullptr;
}

void *operator new[](std::size_t size, const std::nothrow_t &) noexcept {
    return operator new(size, std::nothrow);
}

void operator delete(void *ptr) throw() {
    if (ptr != nullptr) scalable_free(ptr);
}

void operator delete(void *ptr, std::size_t) noexcept {
    if (ptr != nullptr) scalable_free(ptr);
}


void operator delete[](void *ptr) noexcept {
    operator delete(ptr);
}

void operator delete[](void *ptr, std::size_t) noexcept {
    operator delete(ptr);
}
void operator delete(void *ptr, const std::nothrow_t &) noexcept {
    if (ptr != nullptr) scalable_free(ptr);
}

void operator delete[](void *ptr, const std::nothrow_t &) noexcept {
    operator delete(ptr, std::nothrow);
}

