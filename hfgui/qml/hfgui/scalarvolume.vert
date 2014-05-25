#version 330 core
uniform mat4 qt_ModelViewProjectionMatrix;
in vec4 qt_Vertex;
out vec4 entryPoint; // = EntryPoint
out vec4 entryPointTexCoord; // = EntryPoint

void main(void)
{
    gl_Position = qt_ModelViewProjectionMatrix * qt_Vertex;
    entryPoint = qt_Vertex;
    entryPointTexCoord = qt_Vertex; // Should be set to a texture coordinate
}
