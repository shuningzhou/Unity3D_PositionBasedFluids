Shader "PBF/SphereDepth" {
	Properties{
		_mainTex("Texture", 2D) = "white" {}
		_ptlRadius("particle radius", Range(0.1, 1)) = 0.1
	}
	SubShader{
		Pass{

			ZTest LEqual
			Fog{ Mode off }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma enable_d3d11_debug_symbols
			#pragma target 5.0

			#include "Lighting.cginc"
			#include "UnityCG.cginc"

			#pragma geometry geo

			struct Particle
			{
				float4 position;
				float4 velocity;
				float4 predictedPos;
				uint id;
				uint cellId;
				float pad;
				float pad2;
			};

			StructuredBuffer<Particle> _ParticlesBuffer;


			float _ptlRadius;
			Texture2D mainTex;
			
			struct v2g
			{
				float4	pos		: POSITION;
			};

			struct g2f
			{
				float4	pos		: POSITION;
				float2  tex0	: TEXCOORD0;
				float3  eyeSpacePos	: TEXCOORD1;
			};

			struct fout
			{
				float  fragColor	: COLOR;
				float  fragDepth	: DEPTH;
			};

			v2g vert(uint vertex_id : SV_VertexID)
			{
				v2g o;

				//Put this vertex into World space from it's local space
				float4 vPos = _ParticlesBuffer[vertex_id].position;
				o.pos = mul(unity_ObjectToWorld, vPos);
				
				return o;
			}

			[maxvertexcount(4)]
			void geo(point v2g p[1], inout TriangleStream<g2f> triStream)
			{
				float3 look = _WorldSpaceCameraPos - p[0].pos;
				look = normalize(look);

				//use camera's matrix
				float3 up = UNITY_MATRIX_IT_MV[1].xyz;

				float3 right = cross(up, look);

				float halfS = 0.5f * _ptlRadius;

				float4 v[4];
				v[0] = float4(p[0].pos + halfS * right - halfS * up, 1.0f);
				v[1] = float4(p[0].pos + halfS * right + halfS * up, 1.0f);
				v[2] = float4(p[0].pos - halfS * right - halfS * up, 1.0f);
				v[3] = float4(p[0].pos - halfS * right + halfS * up, 1.0f);


				float4x4 vp = UNITY_MATRIX_VP; //put the Quad into View Space
				g2f pIn;
				pIn.pos = mul(vp, v[0]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 0.0f));
				pIn.eyeSpacePos = (mul(UNITY_MATRIX_V, v[0])).xyz;
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[1]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 1.0f));
				pIn.eyeSpacePos = (mul(UNITY_MATRIX_V, v[1])).xyz;
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[2]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 0.0f));
				pIn.eyeSpacePos = (mul(UNITY_MATRIX_V, v[2])).xyz;
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[3]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 1.0f));
				pIn.eyeSpacePos = (mul(UNITY_MATRIX_V, v[3])).xyz;
				triStream.Append(pIn);
			}

			fout frag(g2f input)
			{
				fout OUT;

				float3 N;
				N.xy = input.tex0 * 2.0 - 1.0;
				float r2 = dot(N.xy, N.xy);
				if (r2 > 1.0) discard;	//kill pixels outside circle
				N.z = sqrt(1.0 - r2);

				// calculate depth
				float4 pixelPos = float4 (input.eyeSpacePos + N * _ptlRadius, 1.0);
				float4 clipSpacePos = mul(UNITY_MATRIX_P, pixelPos);

				
				//------------------------------------------------------
				float sphereDepth = clipSpacePos.z / clipSpacePos.w;
				OUT.fragDepth = sphereDepth;

				float color;
				#if defined(UNITY_REVERSED_Z)
				color = 1.0 - sphereDepth;
				#else
				color = sphereDepth;
				#endif

				//if (color > 0.9)
					//color = 0;

				OUT.fragColor = color;

				//OUT.fragColor = sphereDepth;
				//-------------------------------------------------------
				
				//normal test


				/*
				OUT.fragDepth = clipSpacePos.z / clipSpacePos.w;

				float diffuseMul = max(0.0, dot(N, input.lightDir) * 0.5 + 0.5);
				float3 diffuse = input.color * _LightColor0.rgb * diffuseMul;
				float3 ambient = input.color * UNITY_LIGHTMODEL_AMBIENT.rgb;

				OUT.fragColor = float4(ambient + diffuse, 1.0);
				*/

				return OUT;
			}
			ENDCG
		}
	}
}
