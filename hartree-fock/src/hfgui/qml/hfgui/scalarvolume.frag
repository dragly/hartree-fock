//#extension GL_ARB_gpu_shader5 : enable
uniform sampler3D myTexture3D;
uniform highp vec4 ve_eyePosition; // last column of the inverse modelViewMatrix
uniform highp float multiplier;
uniform highp float quality;
uniform highp bool useSquareRootDensity;

varying highp vec4 entryPoint; // = EntryPoint
varying highp vec4 entryPointTexCoord; // = EntryPoint

void main(void)
{
    float stepSize = 1.0 / quality;
    vec3 eye = ve_eyePosition.xyz;
    vec3 exitPoint = eye;
    vec3 direction = exitPoint - entryPoint.xyz;
    direction *= -1;
    float directionLength = length(direction);
    vec3 deltaDir = normalize(direction) * stepSize;
    float deltaDirLength = length(deltaDir);
    vec3 startPoint = entryPointTexCoord.xyz + vec3(0.5, 0.5, 0.5); // TODO: Remove this after implementing proper texture coordinates
    vec3 voxelCoord = startPoint;
    float colorAcummulated = 0.0;
    float redAcummulated = 0.0;
    for(int i = 0; i < int(1.732 / stepSize); i++) { // 1.732 = cube diagonal
        voxelCoord += deltaDir;
        vec4 voxelValue = multiplier * texture3D(myTexture3D, voxelCoord);
        if(useSquareRootDensity) {
            voxelValue.x = sqrt(voxelValue.x);
        }
        colorAcummulated += voxelValue.x * stepSize;
        if(voxelValue.x < 0.3) {
            redAcummulated += voxelValue.x * voxelValue.x * stepSize * 50.0;
        }
        if(voxelCoord.x > 1.0 || voxelCoord.y > 1.0 || voxelCoord.z > 1.0
                || voxelCoord.x < 0.0 || voxelCoord.y < 0.0 || voxelCoord.z < 0.0) {
            break;
        }
    }
    colorAcummulated = clamp(colorAcummulated, 0.0, 1.0);
//    redAcummulated = clamp(redAcummulated, 0.0, 1.0);
    gl_FragColor = vec4(1.0, 1.0, 1.0, colorAcummulated);

}
