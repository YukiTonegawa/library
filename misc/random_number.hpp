#ifndef _RANDOM_NUMBER_H_
#define _RANDOM_NUMBER_H_
#include <random>
#include <cassert>
unsigned long long random_once(){
  static std::random_device seed_gen;
  static std::mt19937_64 engine(seed_gen());
  static unsigned long long ret = engine();
  return ret;
}

unsigned long long random_number(){
  static std::random_device seed_gen;
  static std::mt19937_64 engine(seed_gen());
  return engine();
}

// [low, high]
unsigned long long random_number(unsigned long long low, unsigned long long high){
  static std::random_device seed_gen;
  static std::mt19937_64 engine(seed_gen());
  assert(high >= low);
  unsigned long long diff = high - low + 1;
  return engine() % diff + low;
}
#endif