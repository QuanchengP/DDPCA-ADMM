#ifndef _OctreeElement_hpp
#define _OctreeElement_hpp

#include "../General/General.hpp"

#include <array>
#include <vector>

namespace Ddpca {

/****************************************************************************************************/
//all lines of a hexahedron
inline const std::array<std::array<I64, 2>, 12> hexaLine = {
	{{0, 1}, {1, 2}, {2, 3}, {3, 0}, 
	{0, 4}, {1, 5}, {2, 6}, {3, 7}, 
	{4, 5}, {5, 6}, {6, 7}, {7, 4}}
};

//all faces of a hexahedron, the normal directions to the outside
inline const std::array<std::array<I64, 4>, 6> hexaFace = {
	{{0, 3, 2, 1}, {4, 5, 6, 7}, 
	{0, 4, 7, 3}, {1, 2, 6, 5}, 
	{0, 1, 5, 4}, {3, 7, 6, 2}}
};

/****************************************************************************************************/
// Octree element class for global/local mesh refinement
class OctreeElement {
public:
    // Private data members for better encapsulation
    I64 parent;                                // Index of parent element (-1 if root)
    std::array<I64, 8> cornerNodes;            // Fixed-size array instead of vector to reduce dynamic memory allocation
    I64 level;                                 // Refinement level
    I64 refinementPattern;                     // Refinement pattern as defined in comments
    std::vector<I64> children;                 // Indices of child elements
    
public:
    // Refinement pattern constants for better type safety and readability
    static constexpr I64 REFINEMENT_FULL = 0;       // Split in all three directions (xi, eta, zeta)
    static constexpr I64 REFINEMENT_XI_ETA = 1;     // Split in xi and eta directions
    static constexpr I64 REFINEMENT_ETA_ZETA = 2;   // Split in eta and zeta directions
    static constexpr I64 REFINEMENT_ZETA_XI = 3;    // Split in zeta and xi directions
    static constexpr I64 REFINEMENT_XI = 4;         // Split only in xi direction
    static constexpr I64 REFINEMENT_ETA = 5;        // Split only in eta direction
    static constexpr I64 REFINEMENT_ZETA = 6;       // Split only in zeta direction
    static constexpr I64 REFINEMENT_NONE = 7;       // Not to be refined
    
    // Default constructor - marked noexcept to improve container operation efficiency
    OctreeElement() noexcept : 
        parent(-1), 
        level(0), 
        refinementPattern(REFINEMENT_NONE) 
    {
        // No need for resize or clear, array initializes automatically
    }
    
    // Copy constructor - marked noexcept
    OctreeElement(const OctreeElement& other) noexcept = default;
    
    // Move constructor - marked noexcept to allow containers to use move instead of copy during resizing
    OctreeElement(OctreeElement&& other) noexcept = default;
    
    // Copy assignment operator - marked noexcept
    OctreeElement& operator=(const OctreeElement& other) noexcept = default;
    
    // Move assignment operator - marked noexcept
    OctreeElement& operator=(OctreeElement&& other) noexcept = default;
    
    // Destructor
    ~OctreeElement() noexcept = default;

}; // class OctreeElement

} // namespace Ddpca

#endif // _OctreeElement_hpp