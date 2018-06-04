// Defines basic types used across all commands in the project. All commands
// shall include this, so we always use same types across all commands.
// All classes are templates and they expect FType or SType to be passed, in
// such case we shall pass type aliases defined here.

using F = float;  // default type for floating point number
using S = short;  // default type for pixel or lor index
                  // (2 bytes to save memory)
using Hit = int;  // default type for holding number of hits
