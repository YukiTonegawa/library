#ifndef _TRAITS_H_
#define _TRAITS_H_
#include <cassert>
#include <numeric>
#include <type_traits>

// 符号なし整数 -> 符号付き整数
template <class T>
using make_signed_t = typename std::make_signed<T>::type;
// 符号付き整数 -> 符号なし整数
template <class T>
using make_unsigned_t = typename std::make_unsigned<T>::type;
// 符号付き32bit整数か
template <class T>
using is_signed_int32 = typename std::conditional<std::is_same<T, int>::value || std::is_same<T, int32_t>::value, std::true_type, std::false_type>::type;
// 符号なし32bit整数か
template <class T>
using is_unsigned_int32 = typename std::conditional<std::is_same<T, unsigned int>::value || std::is_same<T, uint32_t>::value, std::true_type, std::false_type>::type;
// 符号付き64bit整数か
template <class T>
using is_signed_int64 = typename std::conditional<std::is_same<T, long long int>::value || std::is_same<T, int64_t>::value, std::true_type, std::false_type>::type;
// 符号なし64bit整数か
template <class T>
using is_unsigned_int64 = typename std::conditional<std::is_same<T, unsigned long long>::value || std::is_same<T, uint64_t>::value, std::true_type, std::false_type>::type;
// 符号付き128bit整数か
template <class T>
using is_signed_int128 = typename std::conditional<std::is_same<T, __int128_t>::value || std::is_same<T, __int128>::value, std::true_type, std::false_type>::type;
// 符号なし128bit整数か
template <class T>
using is_unsigned_int128 = typename std::conditional<std::is_same<T, __uint128_t>::value || std::is_same<T, unsigned __int128>::value, std::true_type, std::false_type>::type;
// 128bit整数か
template <class T>
using is_int128 = typename std::conditional<std::is_same<T, __int128_t>::value || std::is_same<T, __int128>::value || std::is_same<T, __uint128_t>::value || std::is_same<T, unsigned __int128>::value, std::true_type, std::false_type>::type;
// 32bitまたは64bitの整数か
template<class T>
using is_intle64 = typename std::conditional<is_signed_int32<T>::value || is_signed_int64<T>::value || is_unsigned_int32<T>::value || is_unsigned_int64<T>::value, std::true_type, std::false_type>::type;
// 32bitまたは64bitの符号付き整数か
template<class T>
using is_signed_intle64 = typename std::conditional<is_signed_int32<T>::value || is_signed_int64<T>::value, std::true_type, std::false_type>::type;
// 32bitまたは64bitの符号なし整数か
template<class T>
using is_unsigned_intle64 = typename std::conditional<is_unsigned_int32<T>::value || is_unsigned_int64<T>::value, std::true_type, std::false_type>::type;

template<class T>
using is_not_t = typename std::conditional<T::value, std::false_type, std::true_type>::type;
#endif