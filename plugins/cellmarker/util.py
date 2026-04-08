import os


async def get_query_res(request):
    from starlette.responses import FileResponse, Response

    query_res_name = request.path_params["query_res_name"]
    query_res_path = f"./query_results/{query_res_name}"

    if not os.path.isfile(query_res_path):
        return Response(
            content={"error": "query_res not found"}, media_type="application/json"
        )

    return FileResponse(query_res_path)


def add_query_res_route(server):
    from starlette.routing import Route

    server._additional_http_routes = [
        Route("/query_results/{query_res_name}", endpoint=get_query_res)
    ]
