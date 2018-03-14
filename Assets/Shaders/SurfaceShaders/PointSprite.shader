// Upgrade NOTE: replaced '_Object2World' with 'unity_ObjectToWorld'

// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'
// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'
// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "PBF/PointSprite" {
	Properties {
		_mainTex("Texture", 2D) = "white" {}
		_ptlRadius("particle radius", Range(0.1, 1)) = 0.1
		_color("color", Color) = (0.5, 0.5, 0.5, 1)
	}
	SubShader{
		Pass{

			ZTest LEqual
			Fog{ Mode off }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			//#pragma enable_d3d11_debug_symbols
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
			float4 _color;

			struct v2g
			{
				float4	pos		: POSITION;
				float4  color	: COLOR0;
				float4 lightDir	: TEXCOORD1;
			};

			struct g2f
			{
				float4	pos		: POSITION;
				float2  tex0	: TEXCOORD0;
				float3  eyeSpacePos	: TEXCOORD1;
				float4  color	: COLOR0;
				float4 lightDir	: TEXCOORD2;
			};

			struct fout
			{
				float4 fragColor	: COLOR;
				float fragDepth		: DEPTH;
			};

			v2g vert(uint vertex_id : SV_VertexID)
			{
				v2g o;

				//Put this vertex into World space from it's local space
				float4 vPos = _ParticlesBuffer[vertex_id].position;
				o.pos = mul(unity_ObjectToWorld, vPos);
				o.color = _color;
				o.lightDir = mul(UNITY_MATRIX_MV, float4(normalize(_WorldSpaceLightPos0.xyz), 0.0));

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
				pIn.pos			= mul(vp, v[0]);
				pIn.tex0		= MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 0.0f));
				pIn.color		= p[0].color;
				pIn.eyeSpacePos		= (mul(UNITY_MATRIX_V, v[0])).xyz;
				pIn.lightDir	= p[0].lightDir;
				triStream.Append(pIn);

				pIn.pos			= mul(vp, v[1]);
				pIn.tex0		= MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 1.0f));
				pIn.color		= p[0].color;
				pIn.eyeSpacePos		= (mul(UNITY_MATRIX_V, v[1])).xyz;
				pIn.lightDir	= p[0].lightDir;
				triStream.Append(pIn);

				pIn.pos			= mul(vp, v[2]);
				pIn.tex0		= MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 0.0f));
				pIn.color		= p[0].color;
				pIn.eyeSpacePos		= (mul(UNITY_MATRIX_V, v[2])).xyz;
				pIn.lightDir	= p[0].lightDir;
				triStream.Append(pIn);

				pIn.pos			= mul(vp, v[3]);
				pIn.tex0		= MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 1.0f));
				pIn.color		= p[0].color;
				pIn.eyeSpacePos		= (mul(UNITY_MATRIX_V, v[3])).xyz;
				pIn.lightDir	= p[0].lightDir;
				triStream.Append(pIn);
			}

			//with problems!!!!
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
				OUT.fragDepth = clipSpacePos.z / clipSpacePos.w;

				float diffuseMul = max(0.0, dot(N, input.lightDir) * 0.5 + 0.5);
				float3 diffuse = input.color * _LightColor0.rgb * diffuseMul;
				float3 ambient = input.color * UNITY_LIGHTMODEL_AMBIENT.rgb;

				OUT.fragColor = float4(ambient + diffuse, 1.0);


				//OUT.fragColor = input.color;
				//OUT.fragColor = float4(ambient, 1.0);

				return OUT;

				/*
				float4 lightDir = float4(0.0, 0.0, 1.0, 0.0); //no directed light
				float diffuse = max(0.0, dot(N, lightDir));
				
				return input.color;
				*/
			}


			/*
			struct ps_input {
				float4 pos : SV_POSITION;
			};

			//Our vertex function simply fetches a point from the buffer corresponding to the vertex index
			//which we transform with the view-projection matrix before passing to the pixel program.
			ps_input vert(uint vertex_id : SV_VertexID)
			{
				ps_input o;
				float3 worldPos = _ParticlesBuffer[vertex_id].v3Position;
				o.pos = mul(UNITY_MATRIX_VP, float4(worldPos, 1.0f));
				return o;
			}

			//Pixel function returns a solid color for each point.
			float4 frag(ps_input i) : COLOR
			{
				return float4(1,0.5f,0.0f,1);
			}
			*/



			/*
			struct v2f
			{
				float2 texCoord		: TEXCOORD0;
				float3 eyeSpacePos	: TEXCOORD1;
				float  sphereRadius : TEXCOORD2;
				float4 color		: COLOR0;
				float4 position		: SV_POSITION;
			};

			struct fout
			{
				float4 fragColor	: COLOR;
				float  fragDepth	: DEPTH;
			};

			float _pointScale;

			v2f vert(float4 texcoord : TEXCOORD0, uint instance_id : SV_InstanceID)
			{
				v2f o;
				o.color = float4(normalize(_ParticlesBuffer[instance_id].v3Velocity), 1.0f);
				o.sphereRadius = 0.5f;
				float4 pos = float4(_ParticlesBuffer[instance_id].v3Position, 1.0f);
				o.eyeSpacePos = (mul(UNITY_MATRIX_V, pos)).xyz;
				o.texCoord = texcoord;
				o.position = mul(UNITY_MATRIX_VP, pos);

				return o;
			}


			fout frag(v2f o)
			{
				fout OUT;

				//calculate eye-space sphere normal from texture coordinates
				float3 N;
				N.xy = o.texCoord * 2.0 - 1.0;
				float r2 = dot(N.xy, N.xy);
				if (r2 > 1.0) discard;	//kill pixels outside circle
				N.z = -sqrt(1.0 - r2);

				// calculate depth
				float4 pixelPos = float4 (o.eyeSpacePos + N * o.sphereRadius, 1.0);
				float4 clipSpacePos = mul(pixelPos, UNITY_MATRIX_P);
				OUT.fragDepth = clipSpacePos.z / clipSpacePos.w;

				float4 lightDir = float4(0.0, 0.0, 1.0, 0.0); //no directed light
				float diffuse = max(0.0, dot(N, lightDir));
				OUT.fragColor = diffuse * o.color;

				return OUT;
			}
		*/
			ENDCG
		}
	}
}
