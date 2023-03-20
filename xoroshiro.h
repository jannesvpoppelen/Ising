#ifndef XOROSHIRO_H
#define XOROSHIRO_H

#include <stdint.h>

static inline uint64_t rotl(const uint64_t x, int k);

extern uint64_t s[2];

extern uint64_t next(void);

extern void jump(void);

extern void long_jump(void);

#endif /* XOROSHIRO_H */