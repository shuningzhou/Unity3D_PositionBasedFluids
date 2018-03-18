Shader "PBF/SurfaceShading" {
	Properties {
		_NormalTex("Base (RGB)", 2D) = "white" {}
	}
	SubShader{
		Fog{ Mode off }
		
		GrabPass{ "_GrabTexture" }

		Pass{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma target 5.0
			#pragma enable_d3d11_debug_symbols

			#include "Lighting.cginc"
			#include "UnityCG.cginc"

			sampler2D _NormalTex;
			sampler2D _DepthTex;
			sampler2D _GrabTexture;
			sampler2D _ThicknessTex;

			samplerCUBE _Cube;

			float4 color;
			float  refractionIndex; 
			int shininess;
			float4x4 viewMatrix;
			int isTransparent;

			struct v2f
			{
				float4 pos			: POSITION;
				float2 coord		: TEXCOORD0;
			};

			struct fout
			{
				float4 color	: COLOR;
			};

			v2f vert(appdata_base v)
			{
				v2f o;

				o.pos = UnityObjectToClipPos(v.vertex);
				o.coord = v.texcoord.xy;

				return o;
			}

			float3 uvToEye(float2 texCoord, float z)
			{
				// Convert texture coordinate to view space
				float zfar = _ProjectionParams.z;
				float znear = _ProjectionParams.y;

				float2 unitcubePos = (texCoord * 2.0 - 1.0);
				float a = zfar / (zfar - znear);
				float b = zfar * znear / (znear - zfar);
				float rd = b / (z - a);
				return float3(unitcubePos.x, unitcubePos.y, -1.0) * rd;
			}

			//Schlick's approximation, assume interface1 is air -> refrIdx = 1
			float Fresnel(float3 v, float3 n, float refrIdx)
			{
				float r0 = pow((1 - refrIdx) / (1 + refrIdx), 2.0);
				float cosTheta = dot(n, v);

				return (r0 + (1 - r0) * pow(1 - cosTheta, 5.0)); //reflection coefficient
			}

			//obtain inverse of unity view matrix
			float4x4 UnityIV(float4x4 viewMatrix)
			{
				//assume row major matrix
				float3x3 subIV = transpose(float3x3(viewMatrix[0].xyz, viewMatrix[1].xyz, viewMatrix[2].xyz));
				//translation
				float3 tIV = (-1 * mul(subIV, float3(viewMatrix[0][3], viewMatrix[1][3], viewMatrix[2][3])));
				return (float4x4(float4(subIV[0], tIV.x),
					float4(subIV[1], tIV.y),
					float4(subIV[2], tIV.z),
					float4(0, 0, 0, 1)));
			}

			fout frag(v2f i)
			{
				fout OUT;

				float3 eyeSpaceN = tex2D(_NormalTex, i.coord).xyz;
				float  depth = tex2D(_DepthTex, i.coord);

				//unity built in matrix is not defined for frame buffer objects
				float3 lightDir = normalize(mul(viewMatrix, _WorldSpaceLightPos0)).xyz;
				//float4 lightDir= mul(UNITY_MATRIX_V, _WorldSpaceLightPos0);

				//diffuse
				//float diffuseMul = max(0.0, dot(eyeSpaceN, lightDir));
				float diffuseMul = max(dot(eyeSpaceN, lightDir) * 0.5 + 0.5, 0.0); //wrapped diffuse
				float3 mydiffuse = color.rgb * diffuseMul * _LightColor0.rgb;

				float3 myambient = color.rgb * UNITY_LIGHTMODEL_AMBIENT.rgb;  //ShadeSH9(i.normal);

				//specualr 
				float3 eyeSpacePos = uvToEye(i.coord, depth);
				float3 viewerVector = normalize(-eyeSpacePos); //equals to float3(0,0,0) - myPos, camera position is float3(0,0,0) (since this is view space)

				float3 halfVector = normalize(lightDir + viewerVector);
				float specularMul = pow(max(dot(eyeSpaceN, halfVector), 0.0), shininess);
				float3 myspecular = _LightColor0.rgb * specularMul;

				//cubemap reflection
				float3 reflDir = reflect(-viewerVector, eyeSpaceN);
				//text cube only works correct with world coordinates

				//reflection rate
				float Rr = Fresnel(viewerVector, eyeSpaceN, refractionIndex);
				float3 myreflection = texCUBE(_Cube, mul(UnityIV(viewMatrix), reflDir)) * Rr;

				//color absorption

				float thickness = tex2D(_ThicknessTex, i.coord);

				//8.0, 2.0, 1.0 -> constant coefficient for each color channel
				float3 beerLambert = float3(exp(-8.0*thickness), exp(-2.0*thickness), exp(-1.2*thickness));

				//background
				float2 uv = i.coord;
				float2 uv2 = uv + (eyeSpaceN.xy * thickness * 0.5); //0.5 -> refraction fall off
				//flip coordinates if unity uv is inverse
#if UNITY_UV_STARTS_AT_TOP
				uv.y = 1 - uv.y;
				uv2.y = 1 - uv2.y;
#endif

				float3 bg = tex2D(_GrabTexture, uv).xyz;

				if (depth > 0.999)
				{
					//output background pixel if not part of fluids
					OUT.color = float4(bg, 1.0);
				}
				else
				{
					float3 opaque = (mydiffuse + myambient) * beerLambert + myreflection + myspecular;

					if (isTransparent)
					{
						float tScale = clamp(thickness * 0.5 + 0.5, 0.0, 1.0);
						float refraction = tex2D(_GrabTexture, uv2) * (1 - Rr);
						
						OUT.color = float4(opaque * tScale + refraction * (1 - tScale), 1.0);
					}
					else
					{
						OUT.color = float4(opaque, 1.0);
					}


					//float3 testColor = float3(1, 1, 1);
					//float4 opaque = float4((mydiffuse + myambient) * beerLambert + myreflection + myspecular, 1);
					//float4 opaque = float4((mydiffuse + myambient) * beerLambert + myspecular + myreflection, 1);
					//OUT.color = opaque;

					//refracted background
					//OUT.color = tex2D(_GrabTexture, uv2);//uv + eyeSpaceN.xy * thickness);

					//color absorption test
					//float3 testColor = float3(1, 1, 1);
					//OUT.color = float4(testColor * beerLambert, 1.0);

					//thickness
					//OUT.color = float4(thickness * 0.5 + 0.5, thickness * 0.5 + 0.5, thickness * 0.5 + 0.5, 1);

					//transparancy test
					//OUT.color = float4(opaque.xyz * thickness + tex2D(_GrabTexture, uv2) * (1- thickness), 1.0);

					//pure reflection
					//OUT.color = float4(mydiffuse + myspecular, 1);
				}

				return OUT;
				
			}
			ENDCG
		}	
	}
}
