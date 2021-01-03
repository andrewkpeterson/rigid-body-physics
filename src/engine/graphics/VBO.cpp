#include "VBO.h"

#include "VBOAttribMarker.h"

#include <iostream>

// This will count up the total size of each vertex, based on the maximum offset + numElements
unsigned int calculateFloatsPerVertex(const std::vector<VBOAttribMarker> &markers) {
    unsigned int max = 0;
    for (auto it = markers.begin(); it!= markers.end(); it++) {
        max = std::max(max, static_cast<unsigned int>(it->numElements + it->offset / sizeof(GLfloat)));
    }
    return max;
}

VBO::VBO(const float *data, int sizeInFloats, std::vector<VBOAttribMarker> markers, GEOMETRY_LAYOUT layout) :
    m_handle(-1),
    m_markers(markers),
    m_bufferSizeInFloats(sizeInFloats),
    m_numberOfFloatsOrIntsPerVertex(calculateFloatsPerVertex(markers)),
    m_stride(m_numberOfFloatsOrIntsPerVertex * sizeof(GLfloat)),
    m_triangleLayout(layout)
{
    // Generate vbo
    glGenBuffers(1, &m_handle);

    // Bind and send data to vbo
    glBindBuffer(GL_ARRAY_BUFFER, m_handle);
    glBufferData(GL_ARRAY_BUFFER, sizeInFloats * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

VBO::VBO(const int *data, int sizeInInts, std::vector<VBOAttribMarker> markers, GEOMETRY_LAYOUT layout) :
    m_handle(-1),
    m_markers(markers),
    m_bufferSizeInFloats(sizeInInts),
    m_numberOfFloatsOrIntsPerVertex(calculateFloatsPerVertex(markers)),
    m_stride(m_numberOfFloatsOrIntsPerVertex * sizeof(GLfloat)),
    m_triangleLayout(layout)
{
    // Generate vbo
    glGenBuffers(1, &m_handle);

    // Bind and send data to vbo
    glBindBuffer(GL_ARRAY_BUFFER, m_handle);
    glBufferData(GL_ARRAY_BUFFER, sizeInInts * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// This is called a copy constructor, you don't have to worry about it
VBO::VBO(VBO &&that) :
    m_handle(that.m_handle),
    m_markers(std::move(that.m_markers)),
    m_bufferSizeInFloats(that.m_bufferSizeInFloats),
    m_numberOfFloatsOrIntsPerVertex(that.m_numberOfFloatsOrIntsPerVertex),
    m_stride(that.m_stride),
    m_triangleLayout(that.m_triangleLayout)
{
    that.m_handle = 0;
}

// This will also deep copy a VBO
VBO& VBO::operator=(VBO &&that) {
    this->~VBO();

    m_handle = that.m_handle;
    m_markers = std::move(that.m_markers);
    m_bufferSizeInFloats = that.m_bufferSizeInFloats;
    m_numberOfFloatsOrIntsPerVertex = that.m_numberOfFloatsOrIntsPerVertex;
    m_stride = that.m_stride;
    m_triangleLayout = that.m_triangleLayout;

    that.m_handle = 0;

    return *this;
}

VBO::~VBO()
{
    glDeleteBuffers(1, &m_handle);
}

void VBO::bind() const {
    glBindBuffer(GL_ARRAY_BUFFER, m_handle);
}

void VBO::bindAndEnable() const {
    // Enable vertex attributes, and specify data layout
    bind();
    for (unsigned int i = 0; i < m_markers.size(); i++) {
        VBOAttribMarker am = m_markers[i];

        glEnableVertexAttribArray(am.index);
        glVertexAttribPointer(am.index, am.numElements, am.dataType, am.dataNormalize, m_stride, reinterpret_cast<GLvoid*>(am.offset));
    }
}

void VBO::unbind() const {
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

int VBO::numberOfFloatsPerVertex() const {
    return m_numberOfFloatsOrIntsPerVertex;
}

int VBO::numberOfVertices() const {
    return m_bufferSizeInFloats / m_numberOfFloatsOrIntsPerVertex;
}

VBO::GEOMETRY_LAYOUT VBO::triangleLayout() const {
    return m_triangleLayout;
}
