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
//    direction *= -1;
    float directionLength = length(direction);
    vec3 deltaDir = normalize(direction) * stepSize;
    float deltaDirLength = length(deltaDir);
    vec3 startPoint = entryPointTexCoord.xyz + vec3(0.5, 0.5, 0.5); // TODO: Remove this after implementing proper texture coordinates
    vec3 voxelCoord = startPoint;
    float colorAcummulated = 0.0;
    vec4 currentColor = vec4(1.0, 1.0, 1.0, 0.0);
    vec4 standardColor = vec4(1.0, 1.0, 1.0, 1.0);
    vec4 highlightColor = vec4(0.2, 0.8, 1.0, 1.0);
    float redAcummulated = 0.0;
    for(int i = 0; i < int(1.732 / stepSize); i++) { // 1.732 = cube diagonal
        voxelCoord += deltaDir;
        vec4 voxelValue = multiplier * stepSize * texture3D(myTexture3D, voxelCoord);
        if(useSquareRootDensity) {
            voxelValue.x = sqrt(voxelValue.x);
        }
        colorAcummulated += voxelValue.x * stepSize;
//        float voxelValueClamped = clamp(voxelValue.x, 0.0, 1.0);
        float voxelValueClamped = clamp(voxelValue.x, 0.0, 1.0);
        float mixerValue = clamp(voxelValue.x * 10, 0.0, 1.0);
        vec4 color = standardColor * mixerValue + highlightColor * (1 - mixerValue);
//        vec4 color = standardColor;
        currentColor.a = color.a * voxelValueClamped + currentColor.a * (1 - color.a * voxelValueClamped);
        currentColor.rgb = (color.rgb * color.a * voxelValueClamped
                + currentColor.rgb * currentColor.a * (1 - voxelValueClamped * color.a)) / currentColor.a;
        currentColor = clamp(currentColor, 0.0, 1.0);
//        currentColor.r = clamp(currentColor.r, 0.0, 1.0);
//        currentColor.g = clamp(currentColor.g, 0.0, 1.0);
//        currentColor.b = clamp(currentColor.b, 0.0, 1.0);
//        currentColor.a = clamp(currentColor.a, 0.0, 1.0);

        if(voxelCoord.x > 1.0 || voxelCoord.y > 1.0 || voxelCoord.z > 1.0
                || voxelCoord.x < 0.0 || voxelCoord.y < 0.0 || voxelCoord.z < 0.0) {
            break;
        }
    }
    colorAcummulated = clamp(colorAcummulated, 0.0, 1.0);
//    redAcummulated = clamp(redAcummulated, 0.0, 1.0);
//    gl_FragColor = vec4(1.0, 1.0, 1.0, colorAcummulated);
    currentColor = clamp(currentColor, 0.0, 1.0);
    gl_FragColor = currentColor;
//    gl_FragColor = vec4(1.0, 1.0, 1.0, currentColor.a);

}
