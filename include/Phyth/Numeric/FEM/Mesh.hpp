#ifndef PHYTH_MESH_HPP
#define PHYTH_MESH_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Vector3.hpp"
#include "Node.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace Phyth {
    class Mesh {
    public:
        std::size_t AddNode(const Vector3<Quantity<Meter> > &position) {
            std::size_t idx = nodes_.size();
            nodes_.emplace_back(position, idx);
            return idx;
        }

        std::size_t AddElement(const std::array<std::size_t, 4> &node_indices) {
            elements_.push_back(node_indices);
            return elements_.size() - 1;
        }

        [[nodiscard]] const std::vector<Node> &GetNodes() const { return nodes_; }
        [[nodiscard]] const std::vector<std::array<std::size_t, 4> > &GetElements() const { return elements_; }

        [[nodiscard]] const Node &GetNode(const std::size_t index) const { return nodes_[index]; }
        [[nodiscard]] const std::array<std::size_t, 4> &GetElement(const std::size_t index) const { return elements_[index]; }

        [[nodiscard]] std::size_t GetDOFPerNode() const { return dof_per_node_; }
        [[nodiscard]] std::size_t GetDOFStart(const std::size_t node_index) const { return node_index * dof_per_node_; }
        [[nodiscard]] std::size_t TotalDOF() const { return nodes_.size() * dof_per_node_; }

        static Mesh MakeBox(const Vector3<Quantity<Meter> > &center,
                            const Vector3<Quantity<Meter> > &size,
                            const std::uint32_t nx, const std::uint32_t ny, const std::uint32_t nz) {
            Mesh mesh;

            const auto dx = size.x / nx;
            const auto dy = size.y / ny;
            const auto dz = size.z / nz;

            const auto origin = center - size * 0.5;

            for (std::uint32_t k = 0; k <= nz; ++k)
                for (std::uint32_t j = 0; j <= ny; ++j)
                    for (std::uint32_t i = 0; i <= nx; ++i) {
                        auto pos = origin + Vector3<Quantity<Meter> >(
                                       i * dx,
                                       j * dy,
                                       k * dz
                                   );
                        mesh.AddNode(pos);
                    }

            auto node_index = [&](const std::uint32_t i, const std::uint32_t j, const std::uint32_t k) -> std::size_t {
                return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
            };

            for (std::uint32_t k = 0; k < nz; ++k)
                for (std::uint32_t j = 0; j < ny; ++j)
                    for (std::uint32_t i = 0; i < nx; ++i) {
                        const std::size_t n0 = node_index(i, j, k);
                        const std::size_t n1 = node_index(i + 1, j, k);
                        const std::size_t n2 = node_index(i, j + 1, k);
                        const std::size_t n3 = node_index(i + 1, j + 1, k);
                        const std::size_t n4 = node_index(i, j, k + 1);
                        const std::size_t n5 = node_index(i + 1, j, k + 1);
                        const std::size_t n6 = node_index(i, j + 1, k + 1);
                        const std::size_t n7 = node_index(i + 1, j + 1, k + 1);

                        mesh.AddElement({n0, n1, n3, n7});
                        mesh.AddElement({n0, n3, n2, n7});
                        mesh.AddElement({n0, n2, n6, n7});
                        mesh.AddElement({n0, n6, n4, n7});
                        mesh.AddElement({n0, n4, n5, n7});
                    }

            return mesh;
        }

    private:
        std::vector<Node> nodes_;
        std::vector<std::array<std::size_t, 4> > elements_;
        std::size_t dof_per_node_ = 3;
    };
}

#endif // PHYTH_MESH_HPP
