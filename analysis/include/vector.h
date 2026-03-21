#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include <unordered_map>
#include <array>
#include <vector>

template < typename T >
class Matrix {
public:
    Matrix(int rows, int)
    : data_(rows) {}

    std::unordered_map< int, T >& operator[](int i) {
        return data_[i];
    }

    const std::unordered_map< int, T >& operator[](int i) const {
        return data_[i];
    }

    std::vector< std::unordered_map< int, T > >& data() {
        return data_;
    }

    const std::vector< std::unordered_map< int, T > >& data() const {
        return data_;
    }

private:
    std::vector< std::unordered_map< int, T > > data_;
};

class V3 {
    friend V3 operator*(double, const V3&);

public:
    V3(double _x = 0, double _y = 0, double _z = 0)
    : x(_x), y(_y), z(_z){};

    V3 operator+(const V3& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }
    V3 operator-(const V3& v) const {
        return {x - v.x, y - v.y, z - v.z};
    }

    V3 operator*(double d) const {
        return {d * x, d * y, d * z};
    }

    V3 operator/(double d) const {
        return {x / d, y / d, z / d};
    }

    double dot(const V3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    V3 cross(const V3& v) const {
        return {y * v.z - z * v.y,
                z * v.x - x * v.z,
                x * v.y - y * v.x};
    }

    double length() const {
        return std::sqrt(dot(*this));
    }

    V3 normalize() const {
        double len = length();
        if (len < 1e-20)
            len = 1e-20;
        return (*this) / len;
    }

    double x, y, z;
};

inline V3 operator*(double d, const V3& v) {
    return v * d;
}

inline std::ostream& operator<<(std::ostream& out, V3 const& v) {
    return out << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

#endif /* VECTOR_H */
