#ifndef PHYTH_NODE_HPP
#define PHYTH_NODE_HPP
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Vector3.hpp"

namespace Phyth {
    class Node {
    public:
        Node(const Vector3<Quantity<Meter>> &position, const std::size_t index)
            : position_(position), index_(index) {}

        [[nodiscard]] const Vector3<Quantity<Meter>>& GetPosition() const { return position_; }
        [[nodiscard]] int GetIndex() const { return index_; }

    private:
        Vector3<Quantity<Meter>> position_;
        std::size_t index_;
    };
}

#endif //PHYTH_NODE_HPP
