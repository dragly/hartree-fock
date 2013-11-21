attribute highp vec4 qt_Vertex;
uniform highp mat4 qt_ModelViewProjectionMatrix;

varying highp vec4 entryPoint; // = EntryPoint
varying highp vec4 entryPointTexCoord; // = EntryPoint

void main(void)
{
    gl_Position = qt_ModelViewProjectionMatrix * qt_Vertex;
    entryPoint = qt_Vertex;
    entryPointTexCoord = qt_Vertex; // Should be set to a texture coordinate
    gl_FogFragCoord = gl_FogCoord;
}
